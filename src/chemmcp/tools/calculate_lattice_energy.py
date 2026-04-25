import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.bonding_data import BORN_EXPONENTS, MADELUNG_CONSTANTS

logger = logging.getLogger(__name__)

# Ionic radii in pm (Shannon radii for common coordination numbers)
IONIC_RADII = {
    "Li+": 76, "Na+": 102, "K+": 138, "Rb": 149, "Cs": 167,
    "Be2+": 45, "Mg2+": 72, "Ca2+": 100, "Sr2+": 118, "Ba2+": 135,
    "Al3+": 53.5, "Sc3+": 74.5, "Y3+": 90, "La3+": 103.2,
    "Ti4+": 60.5, "Zr4+": 72, "Hf4+": 71,
    "O2-": 140, "S2-": 184, "F-": 133, "Cl-": 181, "Br-": 196, "I-": 220,
}

# Known lattice energies (experimental) for validation (kJ/mol)
EXPERIMENTAL_LATTICE_ENERGIES = {
    "NaCl": -786, "KCl": -701, "MgO": -3795, "CaO": -3414, "CaF2": -2837,
    "LiF": -1036, "NaBr": -732, "CsCl": -657, "AgCl": -905, "ZnS": -2610,
    "Al2O3": -15916, "BaO": -3029, "SrO": -3217, "BeO": -4443, "CaCl2": -2258,
}


@ChemMCPManager.register_tool
class CalculateLatticeEnergy(BaseTool):
    __version__ = "0.1.0"
    name = "CalculateLatticeEnergy"
    func_name = 'calculate_lattice_energy'
    description = "Calculate or query lattice energy of ionic solids using Born-Landé equation and Kapustinskii approximation."
    implementation_description = "Calculates lattice energy via Born-Landé equation: U = -(N_A · M · z⁺ · z⁻ · e²)/(4πε₀ r₀)(1 - 1/n), where M is Madelung constant, z are ion charges, r₀ is interionic distance, n is Born exponent. Also provides Kapustinskii approximation as fallback and experimental values for comparison."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Lattice Energy", "Solid State Chemistry", "Born-Landé Equation", "Ionic Compounds"]
    required_envs = []

    code_input_sig = [
        ('compound', 'str', 'N/A', 'Ionic compound formula (e.g., NaCl, MgO, CaF2)'),
        ('method', 'str', 'born_landé', 'Calculation method: born_lande, kapustinskii, or both'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Compound formula, e.g., \"NaCl\" or \"MgO\"'),
    ]
    output_sig = [
        ('compound', 'str', 'Compound formula'),
        ('lattice_energy_kj_mol', 'float', 'Calculated lattice energy (kJ/mol, negative = exothermic)'),
        ('method', 'str', 'Calculation method used'),
        ('parameters_used', 'dict', 'Parameters: Madelung constant, charges, interionic distance, Born exponent'),
        ('experimental_value', 'float', 'Experimental/literature value if available'),
        ('description', 'str', 'Explanation of calculation'),
    ]
    
    
    examples = [
        {'code_input': {'compound': 'NaCl', 'method': 'born-lande'}, 'text_input': {'query': 'NaCl'}, 'output': {'compound': 'NaCl', 'lattice_energy_kj_mol': -786, 'description': '...', 'method': '...', 'parameters_used': {...}, 'experimental_value': -787}},
        {'code_input': {'compound': 'MgO', 'method': 'born-lande'}, 'text_input': {'query': 'MgO'}, 'output': {'compound': 'MgO', 'lattice_energy_kj_mol': -3795, 'description': '...', 'method': '...', 'parameters_used': {...}, 'experimental_value': -3795}},
        {'code_input': {'compound': 'CaF2', 'method': 'kapustinskii'}, 'text_input': {'query': 'CaF2'}, 'output': {'compound': 'CaF2', 'lattice_energy_kj_mol': -2838, 'description': '...', 'method': '...', 'parameters_used': {...}, 'experimental_value': -2838}},
    ]
    def _run_base(self, compound: str, method: str = "born_lande") -> dict:
        comp = compound.strip()
        
        # Parse compound into ions
        ion_data = self._parse_compound(comp)
        if ion_data is None:
            available = sorted(EXPERIMENTAL_LATTICE_ENERGIES.keys())
            raise ChemMCPInputError(
                f"Cannot parse compound '{comp}' or data not available. "
                f"Available compounds with full data: {available}. "
                f"You can also provide custom parameters."
            )
        
        cation, anion, z_plus, z_minus, n_cations, n_anions = ion_data
        
        # Get ionic radii
        r_cat = IONIC_RADII.get(cation)
        r_an = IONIC_RADII.get(anion)
        
        if r_cat is None or r_an is None:
            raise ChemMCPInputError(f"Ionic radius data not available for {cation} or {anion}.")
        
        r0 = r_cat + r_an  # interionic distance in pm
        
        # Get Born exponent (average of cation and anion)
        n_cat = BORN_EXPONENTS.get(cation, 9)
        n_an = BORN_EXPONENTS.get(anion, 9)
        n_avg = (n_cat + n_an) / 2
        
        # Get Madelung constant
        structure_key = comp
        M = MADELUNG_CONSTANTS.get(structure_key)
        if M is None:
            # Default based on stoichiometry
            if n_cations == 1 and n_anions == 1:
                M = 1.74756  # NaCl-type default
            else:
                M = 1.75  # approximate
        
        result = {"compound": comp, "ions": {"cation": cation, "anion": anion}, 
                  "charges": {"z_plus": z_plus, "z_minus": z_minus},
                  "interionic_distance_pm": r0}
        
        methods = []
        if method in ("born_lande", "both"):
            U_bl = self._calc_born_lande(z_plus, z_minus, r0, M, n_avg)
            result["born_lande"] = {
                "lattice_energy_kj_mol": U_bl,
                "madelung_constant": M,
                "born_exponent": round(n_avg, 1),
                "formula": f"U = -(Nₐ·M·z⁺·z⁻·e²)/(4πε₀r₀) × (1 - 1/n)",
            }
            methods.append(("Born-Landé", U_bl))
        
        if method in ("kapustinskii", "both"):
            n_ions_total = n_cations + n_anions
            U_kap = self._calc_kapustinskii(n_ions_total, z_plus, z_minus, r0)
            result["kapustinskii"] = {
                "lattice_energy_kj_mol": U_kap,
                "formula": f"U = -(K·ν·|z⁺·z⁻|/r₀)×(1 - d/r₀), K=121.4 kJ·pm/mol",
            }
            methods.append(("Kapustinskii", U_kap))
        
        # Experimental value for comparison
        exp = EXPERIMENTAL_LATTICE_ENERGIES.get(comp)
        if exp is not None:
            result["experimental_kj_mol"] = exp
            if methods:
                calc_val = methods[0][1]
                error_pct = abs((calc_val - exp) / exp * 100)
                result["deviation_from_experiment_percent"] = round(error_pct, 1)
        
        result["method_used"] = method
        result["note"] = (
            f"Lattice energy is always negative (exothermic). "
            f"More negative values indicate stronger ionic bonding and higher melting points."
        )
        return result

    @staticmethod
    def _calc_born_lande(z_plus, z_minus, r0, M, n):
        """Born-Landé equation: U = -(N_A * M * z+ * z- * e^2) / (4*pi*eps_0 * r0) * (1 - 1/n)"""
        N_A = 6.022e23  # Avogadro's number
        e = 1.602e-19   # elementary charge (C)
        eps_0 = 8.854e-12  # vacuum permittivity (F/m)
        r0_m = r0 * 1e-12  # pm -> m
        # U in J/mol, convert to kJ/mol
        U_j_mol = -(N_A * M * z_plus * abs(z_minus) * e**2) / (4 * math.pi * eps_0 * r0_m) * (1 - 1/n)
        return round(U_j_mol / 1000, 1)

    @staticmethod
    def _calc_kapustinskii(n_ions, z_plus, z_minus, r0):
        """Kapustinskii approximation: U = -K * v * |z+ * z-| / r0 * (1 - d/r0), K=121.4 kJ·pm/mol, d=34 pm"""
        K = 121.4  # kJ*pm/mol
        d = 34     # pm
        n_cation_anion_pairs = n_ions / 2  # approximate formula units
        U = -(K * n_cation_anion_pairs * abs(z_plus * z_minus) / r0) * (1 - d / r0)
        return round(U, 1)

    @staticmethod
    def _parse_compound(formula):
        """Parse simple ionic compound formula into ions and charges."""
        # Known compounds database
        known = {
            "NaCl": ("Na+", "Cl-", +1, -1, 1, 1),
            "KCl": ("K+", "Cl-", +1, -1, 1, 1),
            "LiF": ("Li+", "F-", +1, -1, 1, 1),
            "NaBr": ("Na+", "Br-", +1, -1, 1, 1),
            "CsCl": ("Cs+", "Cl-", +1, -1, 1, 1),
            "AgCl": ("Ag+", "Cl-", +1, -1, 1, 1),
            "MgO": ("Mg2+", "O2-", +2, -2, 1, 1),
            "CaO": ("Ca2+", "O2-", +2, -2, 1, 1),
            "SrO": ("Sr2+", "O2-", +2, -2, 1, 1),
            "BaO": ("Ba2+", "O2-", +2, -2, 1, 1),
            "BeO": ("Be2+", "O2-", +2, -2, 1, 1),
            "CaF2": ("Ca2+", "F-", +2, -1, 1, 2),
            "CaCl2": ("Ca2+", "Cl-", +2, -1, 1, 2),
            "MgCl2": ("Mg2+", "Cl-", +2, -1, 1, 2),
            "ZnS": ("Zn2+", "S2-", +2, -2, 1, 1),
            "Al2O3": ("Al3+", "O2-", +3, -2, 2, 3),
        }
        return known.get(formula)


if __name__ == "__main__":
    run_mcp_server()
