import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Physical constants
_R = 8.314462618   # J/(mol·K)

@ChemMCPManager.register_tool
class EntropyCalculatorNew(BaseTool):
    """
    熵计算工具 — 计算第三定律熵、混合熵、反应熵变。
    
    支持不同贡献（平动、转动、振动）的熵分解。
    """
    __version__ = "0.1.0"
    name = "EntropyCalculatorNew"
    func_name = "calculate_entropy"
    description = "Calculate absolute entropy (3rd law), entropy of mixing, reaction entropy change (ΔS), and entropy contributions from translational, rotational, vibrational modes."
    implementation_description = "Uses Sackur-Tetrode for translational entropy, statistical mechanics formulas for rotational/vibrational entropy, ΔS_mixing = -R·Σ(x_i·ln x_i) for ideal mixing, and ΔS_rxn = Σ ν_i·S°(i)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Entropy", "Thermodynamics", "Third Law", "Mixing", "Statistical Mechanics", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("calculation_type", "str", "N/A", "'reaction' (ΔS from standard entropies), 'mixing' (entropy of mixing), 'decomposition' (S_trans + S_rot + S_vib), or 'phase_change'"),
        ("standard_entropies", "list", "N/A", "Standard molar entropies in J/(mol·K). For reaction: [products..., reactants...]. For mixing: [S1, S2, ...]."),
        ("stoich_coeffs", "list", "N/A", "Stoichiometric coefficients (+ products, - reactants). For mixing: mole fractions [x1, x2, ...]."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin."),
        ("molecular_mass_amu", "float", "N/A", "Molecular mass for decomposition mode."),
        ("volume_m3", "float", "0.02445", "Volume in m^3 for decomposition mode."),
        ("rotational_constant_cm", "float", "N/A", "Rotational constant B in cm^-1 for linear molecule decomposition."),
        ("vibrational_frequencies_cm", "list", "[]", "Vibrational frequencies in cm^-1 for decomposition."),
        ("molecule_type", "str", "diatomic", "'diatomic' or 'nonlinear' for decomposition."),
        ("sigma", "int", "1", "Symmetry number."),
        ("delta_H_phase", "float", "N/A", "Enthalpy of phase change J/mol for phase_change mode."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "String: calc_type|[S_values]|[coeffs]|T|[optional params...]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with delta_S or total_S, breakdown by contribution if applicable."),
    ]

    examples = [
        {
            "code_input": {
                "calculation_type": "reaction",
                "standard_entropies": [213.7, 188.8, 205.1, 130.7],
                "stoich_coeffs": [2, 1, -1, -3],
                "temperature_k": 298.15,
            },
            "text_input": {
                "input_str": "reaction|[213.7,188.8,205.1,130.7]|[2,1,-1,-3]|298.15"
            },
            "output": {
                "result": {
                    "delta_S_J_mol_K": -217.6,
                    "description": "N2 + 3H2 → 2NH3",
                }
            },
        },
        {
            "code_input": {
                "calculation_type": "mixing",
                "stoich_coeffs": [0.5, 0.5],
            },
            "text_input": {
                "input_str": "mixing||[0.5,0.5]"
            },
            "output": {
                "result": {
                    "delta_S_mixing_J_mol_K": 5.76,
                    "n_components": 2,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _calc_reaction_entropy(self, S_list: List[float], coeffs: List[float], T: float) -> dict:
        """ΔS_rxn = Σ ν_i · S°_i"""
        if len(S_list) != len(coeffs):
            raise ChemMCPError("Length mismatch between entropies and coefficients.")
        
        dS = sum(c * s for c, s in zip(coeffs, S_list))
        
        return {
            "calculation_type": "reaction",
            "delta_S_J_mol_K": round(dS, 4),
            "delta_S_kJ_mol_K": round(dS / 1000.0, 4),
            "temperature_K": T,
            "order_increase": dS > 0,
        }

    def _calc_mixing_entropy(self, mole_fractions: List[float]) -> dict:
        """ΔS_mix = -R · Σ x_i · ln(x_i)"""
        total = sum(mole_fractions)
        xs = [x / total for x in mole_fractions]
        
        dS_mix = 0.0
        contributions = []
        for i, x in enumerate(xs):
            if x <= 0:
                continue
            term = -x * math.log(x)
            dS_mix += _R * term
            contributions.append({"component": i+1, "mole_fraction": round(x, 4), "contribution_J_mol_K": round(_R * term, 4)})
        
        return {
            "calculation_type": "mixing",
            "delta_S_mixing_J_mol_K": round(dS_mix, 4),
            "n_components": len(xs),
            "mole_fractions": [round(x, 4) for x in xs],
            "contributions": contributions,
        }

    def _calc_decomposition(self, T: float, M_amu: float, V: float,
                             rot_const_cm: float, vib_freqs: list,
                             mol_type: str, sigma: int) -> dict:
        """Break down entropy into translational, rotational, vibrational components."""
        import sys
        _NA = 6.02214076e23
        _KB = 1.380649e-23
        _H = 6.62607015e-34
        _AMU_TO_KG = 1.66054e-27
        
        # Translational entropy (Sackur-Tetrode)
        M_kg = M_amu * _AMU_TO_KG
        Lambda = _H / math.sqrt(2 * math.pi * M_kg * _KB * T)
        q_t = V / Lambda ** 3
        S_trans = _R * (math.log(q_t / _NA) + 2.5)  # per mole
        
        # Rotational entropy
        if mol_type.lower().strip() == "diatomic" and rot_const_cm:
            theta_rot = (_H * 2.99792458e10 * rot_const_cm) / _KB  # B in cm^-1 → theta in K
            if T > 0.1 * theta_rot:
                q_r = T / (sigma * theta_rot)
            else:
                q_r = 1.0
            S_rot = _R * (math.log(max(q_r, 1e-300)) + 1.0)
        else:
            S_rot = 0.0
        
        # Vibrational entropy
        S_vib = 0.0
        vib_contributions = []
        for nu_cm in (vib_freqs or []):
            if nu_cm <= 0:
                continue
            x = (_H * 2.99792458e10 * nu_cm) / (_KB * T)
            if x > 500:
                continue
            ex = math.exp(x)
            s_mode = _R * (x / (ex - 1.0) - math.log(1.0 - math.exp(-x)))
            S_vib += s_mode
            vib_contributions.append({"freq_cm-1": nu_cm, "S_J_mol_K": round(s_mode, 2)})
        
        S_total = S_trans + S_rot + S_vib
        
        return {
            "calculation_type": "decomposition",
            "total_S_J_mol_K": round(S_total, 4),
            "S_trans_J_mol_K": round(S_trans, 4),
            "S_rot_J_mol_K": round(S_rot, 4),
            "S_vib_J_mol_K": round(S_vib, 4),
            "vibrational_contributions": vib_contributions,
            "temperature_K": T,
        }

    def _calc_phase_change(self, dH_phase: float, T: float) -> dict:
        """ΔS_phase = ΔH_phase / T."""
        if T <= 0:
            raise ChemMCPError("Temperature must be positive.")
        dS = dH_phase / T
        
        return {
            "calculation_type": "phase_change",
            "delta_S_J_mol_K": round(dS, 4),
            "delta_H_J_mol": round(dH_phase, 4),
            "temperature_K": T,
        }

    def _run_base(self, calculation_type: str, standard_entropies: List[float] = None,
                  stoich_coeffs: List[float] = None, temperature_k: float = 298.15,
                  molecular_mass_amu: float = None, volume_m3: float = 0.02445,
                  rotational_constant_cm: float = None, vibrational_frequencies_cm: List[float] = None,
                  molecule_type: str = "diatomic", sigma: int = 1,
                  delta_H_phase: float = None) -> dict:
        calc_type = calculation_type.lower().strip()
        
        if calc_type == "reaction":
            if standard_entropies is None or stoich_coeffs is None:
                raise ChemMCPError("'reaction' mode requires standard_entropies and stoich_coeffs.")
            return self._calc_reaction_entropy(standard_entropies, stoich_coeffs, temperature_k)
        elif calc_type == "mixing":
            if stoich_coeffs is None:
                raise ChemMCPError("'mixing' mode requires stoich_coeffs as mole fractions.")
            return self._calc_mixing_entropy(stoich_coeffs)
        elif calc_type == "decomposition":
            if molecular_mass_amu is None:
                raise ChemMCPError("'decomposition' mode requires molecular_mass_amu.")
            return self._calc_decomposition(temperature_k, molecular_mass_amu, volume_m3,
                                            rotational_constant_cm, vibrational_frequencies_cm,
                                            molecule_type, sigma)
        elif calc_type == "phase_change":
            if delta_H_phase is None:
                raise ChemMCPError("'phase_change' mode requires delta_H_phase.")
            return self._calc_phase_change(delta_H_phase, temperature_k)
        else:
            raise ChemMCPError(
                f"Unknown type: '{calculation_type}'. "
                f"Options: 'reaction', 'mixing', 'decomposition', 'phase_change'."
            )

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.split("|")
            calc_type = parts[0].strip()
            import json
            S_vals = json.loads(parts[1]) if len(parts) > 1 and parts[1].strip() else None
            coeffs = json.loads(parts[2]) if len(parts) > 2 and parts[2].strip() else None
            T = float(parts[3]) if len(parts) > 3 else 298.15
            M = float(parts[4]) if len(parts) > 4 and parts[4].strip() else None
            V = float(parts[5]) if len(parts) > 5 else 0.02445
            Br = float(parts[6]) if len(parts) > 6 and parts[6].strip() else None
            Vib = json.loads(parts[7]) if len(parts) > 7 and parts[7].strip() else []
            mt = parts[8] if len(parts) > 8 else "diatomic"
            sig = int(parts[9]) if len(parts) > 9 else 1
            dHp = float(parts[10]) if len(parts) > 10 and parts[10].strip() else None
            return self._run_base(calc_type, S_vals, coeffs, T, M, V, Br, Vib, mt, sig, dHp)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
