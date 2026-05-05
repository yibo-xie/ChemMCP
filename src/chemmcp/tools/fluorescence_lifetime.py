import logging
import math
import json
from typing import Optional, List, Union
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class FluorescenceLifetime(BaseTool):
    """
    荧光寿命（τ）和量子产率（Φ）计算工具。
    
    核心关系式：
      Φ = k_r / (k_r + k_nr) = τ / τ_r
      τ = 1 / (k_r + k_nr)
      τ_r = 1 / k_r  （辐射寿命）
    """
    __version__                 = "0.1.0"
    name                        = "FluorescenceLifetime"
    func_name                   = "calculate_fluorescence"
    description                 = "Calculate fluorescence lifetime (τ) and quantum yield (Φ), and their interrelations with radiative/non-radiative rates."
    implementation_description  = "Uses kinetic relations: τ=1/(kr+knr), Φ=kr/(kr+knr)=τ/τr. Supports forward/inverse calculation from any two known quantities."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Spectroscopy", "Photophysics", "Fluorescence", "Quantum Yield"]
    required_envs               = []

    code_input_sig   = [
        ("tau_ns",                    "float",  "N/A",     "Fluorescence lifetime in nanoseconds (ns). Use None if unknown."),
        ("quantum_yield_phi",         "float",  "N/A",     "Quantum yield Φ (0–1). Use None if unknown."),
        ("radiative_rate_kr_per_s",   "float",  "None",    "Radiative rate constant k_r in s⁻¹. Use None if unknown."),
        ("nonradiative_knr_per_s",    "float",  "None",    "Non-radiative rate constant k_nr in s⁻¹. Use None if unknown."),
    ]

    text_input_sig   = [
        ("input_params",              "str",    "N/A",     "Space-separated: 'tau_ns phi kr knr' (use 'None' for unknowns)."),
    ]

    output_sig       = [
        ("result",                    "dict",   "Dictionary containing tau_ns, quantum_yield_phi, radiative_lifetime_tau_r_ns, kr_1/s, knr_1/s."),
    ]

    examples         = [
        {
            "code_input": {
                "tau_ns":               5.0,
                "quantum_yield_phi":     0.8,
                "radiative_rate_kr_per_s": None,
                "nonradiative_knr_per_s":  None,
            },
            "text_input": {
                "input_params":          "5.0 0.8 None None",
            },
            "output": {
                "result": {
                    "tau_ns": 5.0,
                    "quantum_yield_phi": 0.8,
                    "radiative_lifetime_tau_r_ns": 6.25,
                    "kr_per_s": 160000000.0,
                    "knr_per_s": 40000000.0,
                }
            },
        },
    ]

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        tau_ns: Optional[float],
        quantum_yield_phi: Optional[float],
        radiative_rate_kr_per_s: Optional[float] = None,
        nonradiative_knr_per_s: Optional[float] = None,
    ) -> dict:
        """
        Core logic: calculate fluorescence parameters from any two known values.
        """
        # Helper to resolve "None" string -> Python None
        def resolve(v):
            if v is None or (isinstance(v, str) and v.strip().lower() == "none"):
                return None
            return float(tau_ns)

        # --- Validate at least two inputs provided ---
        known = sum(1 for v in [tau_ns, quantum_yield_phi, radiative_rate_kr_per_s, nonradiative_knr_per_s]
                     if v is not None and not (isinstance(v, str) and v.strip().lower() == "none"))
        if known < 2:
            raise ChemMCPError("At least two of [tau_ns, phi, kr, knr] must be provided.")

        # Convert ns -> s for internal calc
        tau_s = float(tau_ns) * 1e-9 if tau_ns is not None else None
        kr = float(radiative_rate_kr_per_s) if radiative_rate_kr_per_s is not None else None
        knr = float(nonradiative_knr_per_s) if nonradiative_knr_per_s is not None else None
        phi = float(quantum_yield_phi) if quantum_yield_phi is not None else None

        # --- Validate phi range ---
        if phi is not None and not (0 <= phi <= 1):
            raise ChemMCPError("Quantum yield φ must be between 0 and 1.")

        # --- Derive missing quantities ---
        # Case 1: τ and Φ given → derive kr, knr
        if tau_s is not None and phi is not None:
            if tau_s == 0:
                raise ChemMCPError("τ cannot be zero.")
            total_rate = 1.0 / tau_s
            kr = phi * total_rate
            knr = total_rate - kr
        # Case 2: kr and knr given → derive τ, Φ
        elif kr is not None and knr is not None:
            total_rate = kr + knr
            if total_rate == 0:
                raise ChemMCPError("Total rate (kr+knr) cannot be zero.")
            tau_s = 1.0 / total_rate
            phi = kr / total_rate
        # Case 3: τ and kr given → derive knr, Φ
        elif tau_s is not None and kr is not None:
            total_rate = 1.0 / tau_s
            knr = total_rate - kr
            if knr < 0:
                raise ChemMCPError(f"Derived knr={knr} is negative. Check inputs: kr={kr} cannot exceed total rate {total_rate}.")
            phi = kr / total_rate
        # Case 4: τ and knr given → derive kr, Φ
        elif tau_s is not None and knr is not None:
            total_rate = 1.0 / tau_s
            kr = total_rate - knr
            if kr < 0:
                raise ChemMCPError(f"Derived kr={kr} is negative. Check inputs: knr={knr} cannot exceed total rate {total_rate}.")
            phi = kr / total_rate
        # Case 5: Φ and kr given → derive τ, knr
        elif phi is not None and kr is not None:
            if kr == 0:
                raise ChemMCPError("kr cannot be zero when deriving from φ.")
            tau_r_s = 1.0 / kr
            tau_s = phi * tau_r_s
            knr = kr * (1.0 - phi) / phi if phi > 0 else 0.0
        # Case 6: Φ and knr given → derive kr, τ
        elif phi is not None and knr is not None:
            if phi >= 1:
                # Perfect QY means knr should be ~0
                if knr != 0:
                    raise ChemMCPError("φ=1 implies knr=0, but knr was specified as non-zero.")
                kr = 1.0  # arbitrary; τ depends on additional info
                tau_s = 1.0 / kr
            elif phi == 0:
                kr = 0.0
                tau_s = 1.0 / knr if knr > 0 else float('inf')
            else:
                # φ = kr/(kr+knr) => kr = φ*knr/(1-φ)
                kr = phi * knr / (1.0 - phi)
                total_rate = kr + knr
                tau_s = 1.0 / total_rate
        else:
            raise ChemMCPError("Insufficient information to derive all quantities.")

        # --- Compute derived quantities ---
        tau_final_ns = tau_s * 1e9 if tau_s is not None else None
        tau_r_ns = (1.0 / kr) * 1e9 if (kr is not None and kr > 0) else None

        result = {
            "tau_ns":                  round(tau_final_ns, 6) if tau_final_ns is not None else None,
            "quantum_yield_phi":       round(phi, 6) if phi is not None else None,
            "radiative_lifetime_tau_r_ns": round(tau_r_ns, 6) if tau_r_ns is not None else None,
            "kr_per_s":                round(kr, 4) if kr is not None else None,
            "knr_per_s":                round(knr, 4) if knr is not None else None,
        }

        logger.info(f"Fluorescence calculation: τ={tau_final_ns}ns, Φ={phi}, kr={kr}, knr={knr}")
        return result

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.split()
            if len(parts) < 2:
                raise ValueError("Need at least 2 params: 'tau_ns phi [kr] [knr]'")

            def parse_val(s):
                if s.strip().lower() == "none":
                    return None
                return float(s)

            vals = [parse_val(p) for p in parts]
            tau_ns = vals[0]
            phi = vals[1]
            kr = vals[2] if len(vals) > 2 else None
            knr = vals[3] if len(vals) > 3 else None

            return self._run_base(tau_ns, phi, kr, knr)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'tau_ns phi [kr] [knr]', use 'None' for unknowns.")
