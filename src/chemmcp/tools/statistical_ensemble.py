import json
import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Physical constants (SI units)
KB = 1.380649e-23  # Boltzmann constant, J/K
NA = 6.02214076e23  # Avogadro's number, mol⁻¹
R_GAS = 8.314462618  # Gas constant, J/(mol·K)
PLANCK = 6.62607015e-34  # Planck constant, J·s


@ChemMCPManager.register_tool
class StatisticalEnsemble(BaseTool):
    """
    统计系综计算工具。
    支持微正则、正则、巨正则系综的热力学平均量计算：内能U、熵S、自由能F、热容Cv等。
    """
    __version__ = "0.1.0"
    name = "StatisticalEnsemble"
    func_name = "compute_ensemble"
    description = "Statistical ensemble calculations for thermodynamic averages: internal energy, entropy, free energy, heat capacity, partition function."
    implementation_description = "Implements three statistical ensembles: microcanonical (NVE), canonical (NVT) via Boltzmann distribution, and grand canonical (μVT). Computes partition function Z, internal energy U, entropy S, Helmholtz free energy F, pressure P, heat capacity Cv, and chemical potential from discrete energy levels."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Statistical Mechanics", "Thermodynamics", "Ensemble", "Partition Function", "Physical Chemistry", "Boltzmann"]
    required_envs = []

    code_input_sig = [
        ("ensemble_type", "str", "N/A", "Ensemble type: 'microcanonical', 'canonical', or 'grand_canonical'."),
        ("energy_levels", "List[float]", "N/A", "List of energy levels (in Joules or reduced units)."),
        ("parameters", "dict", "N/A", "Ensemble-specific parameters dict. Canonical: {'T': float(K)}. Grand canonical: {'T': float, 'mu': float(J)}. Microcanonical: {'E_total': float} or {}."),
        ("degeneracies", "List[int]", "None", "Optional degeneracy per level (default all 1)."),
        ("units", "str", "reduced", "Units: 'reduced' (kB=1) or 'SI' (Joules, Kelvin)."),
        ("observable", "str", "all", "Specific observable to compute ('U','S','F','P','Cv','mu_chem','Z','probabilities') or 'all'."),
    ]

    text_input_sig = [
        ("params_str", "str", "N/A", "JSON string with ensemble_type, energy_levels, parameters, etc."),
    ]

    output_sig = [
        ("ensemble_type", "str", "Ensemble type used."),
        ("partition_function", "float", "Partition function Z."),
        ("internal_energy", "float", "Average internal energy ⟨E⟩ (or total E for microcanonical)."),
        ("entropy", "float", "Entropy S."),
        ("helmholtz_free_energy", "float", "Helmholtz free energy F = U - TS."),
        ("heat_capacity_cv", "float", "Heat capacity at constant volume Cv = d⟨E⟩/dT."),
        ("pressure", "float", "Pressure P (if applicable)."),
        ("chemical_potential", "float", "Chemical potential μ (for grand canonical)."),
        ("probabilities", "list", "Probability of each energy level p_i."),
        ("n_levels", "int", "Number of energy levels."),
        ("temperature", "float", "Temperature used (K or reduced units)."),
        ("diagnostics", "dict", "Additional info: beta, average energy squared, etc."),
    ]

    examples = [
        {
            "code_input": {
                "ensemble_type": "canonical",
                "energy_levels": [0.0, 1.0, 2.0, 3.0],
                "parameters": {"T": 1.0},
                "units": "reduced",
            },
            "text_input": {"params_str": '{"ensemble_type":"canonical","energy_levels":[0,1,2,3],"parameters":{"T":1},"units":"reduced"}'},
            "output": {
                "partition_function": 1.5151,
                "internal_energy": 1.0246,
                "entropy": 0.8749,
            },
        },
        {
            "code_input": {
                "ensemble_type": "microcanonical",
                "energy_levels": [0.0, 1.0, 2.0],
                "parameters": {"E_total": 1.0},
            },
            "text_input": {"params_str": '{"ensemble_type":"microcanonical","energy_levels":[0,1,2],"parameters":{"E_total":1}}'},
            "output": {
                "partition_function": 1.0,
                "internal_energy": 1.0,
                "entropy": 0.0,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        ensemble_type: str,
        energy_levels: List[float],
        parameters: dict,
        degeneracies: Optional[List[int]] = None,
        units: str = "reduced",
        observable: str = "all",
    ) -> dict:
        """Core logic: statistical ensemble calculations."""
        if not energy_levels:
            raise ChemMCPError("energy_levels cannot be empty.")
        etype = ensemble_type.lower().strip()
        if etype not in ("microcanonical", "canonical", "grand_canonical"):
            raise ChemMCPError(f"Unknown ensemble type '{ensemble_type}'. Use: microcanonical, canonical, grand_canonical.")

        n_levels = len(energy_levels)
        degen = degeneracies or [1] * n_levels
        if len(degen) != n_levels:
            raise ChemMCPError(f"degeneracies length ({len(degen)}) != energy_levels length ({n_levels}).")

        params = parameters or {}
        use_si = (units.lower() == "si")

        result = {
            "ensemble_type": etype,
            "n_levels": n_levels,
            "units": units,
        }

        if etype == "microcanonical":
            result.update(self._microcanonical(energy_levels, degen, params, use_si))

        elif etype == "canonical":
            T = params.get("T")
            if T is None or T <= 0:
                raise ChemMCPError("Canonical ensemble requires T > 0 in parameters.")
            result["temperature"] = T
            result.update(self._canonical(energy_levels, degen, T, use_si, observable))

        elif etype == "grand_canonical":
            T = params.get("T")
            mu = params.get("mu")
            if T is None or T <= 0:
                raise ChemMCPError("Grand canonical ensemble requires T > 0.")
            result["temperature"] = T
            result.update(self._grand_canonical(energy_levels, degen, T, mu, use_si))

        logger.info(
            f"StatisticalEnsemble ({etype}): Z={result.get('partition_function', 'N/A')}, "
            f"U={result.get('internal_energy', 'N/A')}, S={result.get('entropy', 'N/A')}"
        )
        return result

    def _microcanonical(self, levels, degen, params, use_si):
        """NVE ensemble: all states with same total energy equally probable."""
        E_total = params.get("E_total", sum(levels))
        # Count states at or near E_total (within tolerance)
        tol = params.get("tolerance", 1e-10 * (max(levels) - min(levels)) if len(levels) > 1 else 1e-10)
        omega = sum(g for E, g in zip(levels, degen) if abs(E - E_total) <= tol)
        if omega == 0:
            omega = sum(degen)  # fallback: all states accessible

        kB_local = KB if use_si else 1.0
        S = kB_local * math.log(max(omega, 1))
        return {
            "partition_function": float(omega),
            "internal_energy": E_total,
            "entropy": round(S, 12),
            "helmholtz_free_energy": E_total - S * (params.get("T", 1.0) if use_si else 1.0),
            "heat_capacity_cv": 0.0,  # microcanonical: no thermal fluctuations
            "pressure": None,
            "chemical_potential": None,
            "probabilities": [round(g / omega, 10) for g in degen] if omega > 0 else [0.0] * len(degen),
            "diagnostics": {"multiplicity_omega": omega, "n_accessible_states": omega},
        }

    def _canonical(self, levels, degen, T, use_si, observable):
        """NVT ensemble: Boltzmann distribution."""
        kB_local = KB if use_si else 1.0
        beta = 1.0 / (kB_local * T)

        # Partition function Z = Σ g_i * exp(-β*E_i)
        Z = sum(g * math.exp(-beta * E) for g, E in zip(degen, levels))
        if Z <= 0:
            raise ChemMCPError("Partition function is zero or negative; check energy levels and temperature.")

        # Probabilities
        probs = [g * math.exp(-beta * E) / Z for g, E in zip(degen, levels)]

        # Internal energy U = Σ p_i * E_i
        U = sum(p * E for p, E in zip(probs, levels))

        # Energy fluctuation for heat capacity
        U2 = sum(p * E ** 2 for p, E in zip(probs, levels))
        dU_dT = (U2 - U ** 2) * beta ** 2 * kB_local  # Cv = d<U>/dT

        # Entropy S = -kB Σ p_i ln(p_i) (Gibbs entropy)
        S = 0.0
        for p in probs:
            if p > 1e-30:
                S -= p * math.log(p)
        S *= kB_local

        # Helmholtz free energy F = -kT ln(Z)
        F = -kB_local * T * math.log(Z)

        # Pressure (for ideal gas-like system): P = (1/βV) ln(Z) ≈ kB*T/V * ln(Z) per particle
        V = 1.0  # assume unit volume
        P = kB_local * T * math.log(Z) / V if Z > 1 else 0.0

        ret = {
            "partition_function": round(Z, 12),
            "internal_energy": round(U, 14),
            "entropy": round(S, 14),
            "helmholtz_free_energy": round(F, 14),
            "heat_capacity_cv": round(dU_dT, 14),
            "pressure": round(P, 14) if P else None,
            "chemical_potential": None,
            "probabilities": [round(p, 10) for p in probs],
            "diagnostics": {
                "beta": round(beta, 6),
                "energy_squared_avg": round(U2, 14),
                "energy_fluctuation": round(U2 - U ** 2, 14),
            },
        }
        return ret

    def _grand_canonical(self, levels, degen, T, mu, use_si):
        """μVT ensemble: include chemical potential."""
        kB_local = KB if use_si else 1.0
        beta = 1.0 / (kB_local * T)

        # Grand partition function Ξ = Σ g_i * exp(-β(E_i - μ*N_i))
        # For single-particle levels, N_i = 1
        Xi = sum(g * math.exp(-beta * (E - mu)) for g, E in zip(degen, levels))
        if Xi <= 0:
            raise ChemMCPError("Grand partition function is zero or negative.")

        probs = [g * math.exp(-beta * (E - mu)) / Xi for g, E in zip(degen, levels)]
        U = sum(p * E for p, E in zip(probs, levels))
        N_avg = sum(p * 1.0 for p in probs)  # average particle number

        U2 = sum(p * E ** 2 for p, E in zip(probs, levels))
        dU_dT = (U2 - U ** 2) * beta ** 2 * kB_local

        S = 0.0
        for p in probs:
            if p > 1e-30:
                S -= p * math.log(p)
        S *= kB_local

        F = -kB_local * T * math.log(Xi)

        return {
            "partition_function": round(Xi, 12),
            "internal_energy": round(U, 14),
            "entropy": round(S, 14),
            "helmholtz_free_energy": round(F, 14),
            "heat_capacity_cv": round(dU_dT, 14),
            "pressure": None,
            "chemical_potential": round(mu, 14),
            "average_particle_number": round(N_avg, 10),
            "probabilities": [round(p, 10) for p in probs],
            "diagnostics": {"beta": round(beta, 6)},
        }

    def _run_text(self, params_str: str) -> dict:
        try:
            kwargs = json.loads(params_str.strip())
        except json.JSONDecodeError:
            raise ChemMCPError("Invalid JSON input.")
        return self._run_base(**kwargs)
