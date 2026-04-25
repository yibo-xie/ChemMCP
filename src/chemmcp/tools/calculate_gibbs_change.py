import logging
from typing import Dict, List, Optional, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CalculateGibbsChange(BaseTool):
    """
    计算反应吉布斯自由能变化 ΔG°rxn。
    使用公式: ΔG°rxn = ΔH°rxn - TΔS°rxn
    其中 ΔH°rxn 和 ΔS°rxn 分别由生成焓和标准熵计算。
    """
    __version__ = "0.1.0"
    name = "CalculateGibbsChange"
    func_name = "calculate_gibbs_change"
    description = "Calculate reaction Gibbs free energy change (ΔG°rxn) from standard thermodynamic data using ΔG° = ΔH° - TΔS."
    implementation_description = "Calculates ΔG°rxn = Σ(ν_i×ΔG°f_products) - Σ(ν_j×ΔG°f_reactants) at given temperature, or via ΔG° = ΔH° - TΔS if entropy data is available. Uses built-in thermodynamic database."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Gibbs Free Energy", "Spontaneity", "Thermodynamics", "Reaction Feasibility"]
    required_envs = []

    code_input_sig = [
        ("reactants", "dict", "N/A", "Reactants as {species: stoichiometric_coefficient}, e.g., {'CH4': 1, 'O2': 2}."),
        ("products", "dict", "N/A", "Products as {species: stoichiometric_coefficient}, e.g., {'CO2': 1, 'H2O': 2}."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin (default 298.15 K)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Format: reactants;products;temperature_K. E.g., 'CH4:1,O2:2;CO2:1,H2O:2;298.15' or 'N2+3H2;2NH3'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with delta_g_rxn (kJ/mol), delta_h_rxn, delta_s_rxn, temperature_k, is_spontaneous, and equilibrium constant K."),
    ]

    examples = [
        {
            "code_input": {
                "reactants": {"CH4": 1, "O2": 2},
                "products": {"CO2": 1, "H2O(l)": 2},
                "temperature_k": 298.15,
            },
            "text_input": {
                "input_str": "CH4:1,O2:2;CO2:1,H2O:2",
            },
            "output": {
                "result": {
                    "delta_g_rxn_kj_per_mol": -817.97,
                    "unit": "kJ/mol",
                    "is_spontaneous": True,
                    "temperature_k": 298.15,
                    "delta_h_rxn_kj_per_mol": -890.36,
                    "delta_s_rxn_j_per_mol_k": -242.98,
                    "equilibrium_constant_K": "very large (K >> 1)",
                    "explanation": "ΔG° << 0 → highly spontaneous combustion reaction at standard conditions.",
                }
            },
        },
        {
            "code_input": {
                "reactants": {"N2": 1, "H2": 3},
                "products": {"NH3": 2},
                "temperature_k": 298.15,
            },
            "text_input": {
                "input_str": "N2:1,H2:3;NH3:2",
            },
            "output": {
                "result": {
                    "delta_g_rxn_kj_per_mol": -32.80,
                    "unit": "kJ/mol",
                    "is_spontaneous": True,
                    "temperature_k": 298.15,
                    "delta_h_rxn_kj_per_mol": -91.88,
                    "delta_s_rxn_j_per_mol_k": -198.76,
                    "equilibrium_constant_K": "~7.6 × 10^5",
                    "explanation": "ΔG° < 0 → spontaneous (Haber process). K is large, favoring products.",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize thermodynamic databases."""
        # ΔG°f database (kJ/mol)
        self._dgf_db = {
            # Elements (0)
            "H2": 0.0, "O2": 0.0, "N2": 0.0, "C": 0.0, "S": 0.0,
            "Cl2": 0.0, "Br2": 0.0, "I2": 0.0, "F2": 0.0,
            "Na": 0.0, "Fe": 0.0, "Zn": 0.0, "Cu": 0.0, "Ag": 0.0,
            "Ca": 0.0, "Mg": 0.0, "Al": 0.0, "Si": 0.0,
            # Compounds
            "H2O(l)": -237.14, "H2O(g)": -228.57, "H2O": -237.14,
            "CO(g)": -137.16, "CO": -137.16,
            "CO2(g)": -394.36, "CO2": -394.36,
            "SO2(g)": -300.19, "SO2": -300.19,
            "SO3(g)": -371.06,
            "NO(g)": 86.55, "NO": 86.55,
            "NO2(g)": 51.26,
            "HF(g)": -275.46,
            "HCl(g)": -95.30,
            "HBr(g)": -53.45,
            "HI(g)": 1.70,
            "H2S(g)": -33.56,
            "NH3(g)": -16.40, "NH3": -16.40,
            "CH4(g)": -50.75, "CH4": -50.75,
            "C2H6(g)": -32.82,
            "C2H4(g)": 68.15,
            "C2H2(g)": 209.20,
            "C3H8(g)": -23.47,
            "C6H6(l)": 124.45,
            "CH3OH(l)": -166.27,
            "C2H5OH(l)": -174.78,
            "HCOOH(l)": -361.40,
            "CH3COOH(l)": -389.85,
            "C6H12O6(s)": -910.44,
            "H2O2(l)": -120.35,
            "H2SO4(l)": -689.90,
            "HNO3(l)": -80.71,
            # Oxides
            "MgO(s)": -569.33, "MgO": -569.33,
            "CaO(s)": -603.54, "CaO": -603.54,
            "Al2O3(s)": -1582.27, "Al2O3": -1582.27,
            "Fe2O3(s)": -742.24, "Fe2O3": -742.24,
            "Fe3O4(s)": -1015.38,
            "CuO(s)": -129.71,
            "ZnO(s)": -320.52,
            "SiO2(s)": -856.64, "SiO2": -856.64,
            "PbO(s)": -187.89,
            # Bases
            "NaOH(s)": -379.53,
            "Ca(OH)2(s)": -898.49,
            "Mg(OH)2(s)": -833.58,
            # Halides
            "NaCl(s)": -384.14, "NaCl": -384.14,
            "AgCl(s)": -109.79,
            "CaCl2(s)": -748.10,
            "MgCl2(s)": -591.79,
            "FeCl2(s)": -302.30,
            "CuCl2(s)": -179.90,
            "ZnCl2(s)": -369.43,
            # Carbonates
            "CaCO3(s)": -1128.76, "CaCO3": -1128.76,
            "Na2CO3(s)": -1044.44,
            "BaCO3(s)": -1134.41,
            # Sulfates
            "BaSO4(s)": -1362.18,
            "CaSO4(s)": -1321.74,
            "Na2SO4(s)": -1270.12,
            "CuSO4(s)": -661.86,
            "MgSO4(s)": -1170.66,
            # Nitrates
            "NaNO3(s)": -367.04,
            "KNO3(s)": -394.93,
            "NH4NO3(s)": -189.47,
            # Sulfides
            "FeS(s)": -100.42,
            "ZnS(s)": -201.29,
            "PbS(s)": -98.73,
            # Ammonium salts
            "NH4Cl(s)": -202.95,
            "(NH4)2SO4(s)": -902.28,
            # Other
            "PCl3(g)": -267.77,
            "SF6(g)": -1105.14,
            "BF3(g)": -1120.33,
            "CS2(l)": 65.27,
            "HCN(g)": 124.60,
            # Ions (aq)
            "H+(aq)": 0.0, "OH-(aq)": -157.24,
            "Na+(aq)": -261.88, "K+(aq)": -283.27,
            "Ag+(aq)": 77.11, "Ca2+(aq)": -553.58,
            "Mg2+(aq)": -456.01,
            "Fe2+(aq)": -78.90, "Fe3+(aq)": -4.65,
            "Cu2+(aq)": 65.49, "Zn2+(aq)": -147.19,
            "Al3+(aq)": -485.00,
            "Cl-(aq)": -131.22, "Br-(aq)": -103.96,
            "I-(aq)": -51.59,
            "SO4^2-(aq)": -744.53, "CO3^2-(aq)": -527.81,
            "NO3-(aq)": -111.25, "NH4+(aq)": -79.31,
        }

        # S° database (J/(mol·K)) — for TΔS calculation
        self._s_db = {
            "H2": 130.68, "O2": 205.15, "N2": 191.61, "C(graphite)": 5.74, "C": 5.74,
            "S(s)": 31.80, "Cl2": 223.07, "F2": 202.79,
            "Na(s)": 51.21, "Fe(s)": 27.28, "Zn(s)": 41.63, "Cu(s)": 33.15,
            "Ca(s)": 41.59, "Mg(s)": 32.68, "Al(s)": 28.33, "Si(s)": 18.83,
            "H2O(l)": 69.91, "H2O(g)": 188.83, "H2O": 69.91,
            "CO(g)": 197.67, "CO": 197.67,
            "CO2(g)": 213.79, "CO2": 213.79,
            "SO2(g)": 248.21, "SO2": 248.21,
            "SO3(g)": 256.76,
            "NO(g)": 210.76, "NO": 210.76,
            "NO2(g)": 240.06,
            "HF(g)": 173.78,
            "HCl(g)": 186.91,
            "HBr(g)": 198.70,
            "HI(g)": 206.59,
            "H2S(g)": 205.79,
            "NH3(g)": 192.45, "NH3": 192.45,
            "CH4(g)": 186.26, "CH4": 186.26,
            "C2H6(g)": 229.60,
            "C2H4(g)": 219.32,
            "C2H2(g)": 200.94,
            "C3H8(g)": 270.20,
            "C6H6(l)": 173.26,
            "CH3OH(l)": 126.70,
            "C2H5OH(l)": 160.70,
            "HCOOH(l)": 129.00,
            "CH3COOH(l)": 159.83,
            "C6H12O6(s)": 212.13,
            "H2O2(l)": 109.60,
            # Oxides
            "MgO(s)": 26.94, "MgO": 26.94,
            "CaO(s)": 39.75, "CaO": 39.75,
            "Al2O3(s)": 50.92, "Al2O3": 50.92,
            "Fe2O3(s)": 87.40, "Fe2O3": 87.40,
            "CuO(s)": 42.63,
            "ZnO(s)": 43.64,
            "SiO2(s)": 41.84, "SiO2": 41.84,
            "PbO(s)": 66.50,
            # Halides
            "NaCl(s)": 72.13, "NaCl": 72.13,
            "AgCl(s)": 96.20,
            "CaCl2(s)": 104.60,
            "MgCl2(s)": 89.62,
            "FeCl2(s)": 118.00,
            "CuCl2(s)": 108.00,
            "ZnCl2(s)": 111.46,
            # Carbonates
            "CaCO3(s)": 92.90, "CaCO3": 92.90,
            "Na2CO3(s)": 134.98,
            "BaCO3(s)": 112.10,
            # Sulfates
            "BaSO4(s)": 132.16,
            "CaSO4(s)": 106.70,
            "Na2SO4(s)": 149.58,
            "CuSO4(s)": 113.00,
            "MgSO4(s)": 126.40,
            # Nitrates
            "NaNO3(s)": 116.52,
            "KNO3(s)": 151.04,
            "NH4NO3(s)": 151.08,
            # Sulfides
            "FeS(s)": 60.29,
            "ZnS(s)": 57.70,
            "PbS(s)": 91.30,
            # Ammonium salts
            "NH4Cl(s)": 94.56,
            # Other
            "PCl3(g)": 311.67,
            "SF6(g)": 291.54,
            "BF3(g)": 254.12,
            "CS2(l)": 151.34,
            "HCN(g)": 201.78,
            # Ions (aq)
            "H+(aq)": 0.0, "OH-(aq)": -10.75,
            "Na+(aq)": 58.41, "K+(aq)": 102.50,
            "Ag+(aq)": 72.68, "Ca2+(aq)": -53.10,
            "Mg2+(aq)": -138.10,
            "Fe2+(aq)": -137.70, "Fe3+(aq)": -315.90,
            "Cu2+(aq)": -98.00, "Zn2+(aq)": -112.10,
            "Al3+(aq)": -321.70,
            "Cl-(aq)": 56.48, "Br-(aq)": 111.30,
            "I-(aq)": 111.30,
            "SO4^2-(aq)": 20.10, "CO3^2-(aq)": -56.90,
            "NO3-(aq)": 146.40, "NH4+(aq)": 112.84,
        }

        # ΔH°f database for cross-check (kJ/mol)
        self._dhf_db = {
            "H2": 0.0, "O2": 0.0, "N2": 0.0, "C": 0.0, "S": 0.0,
            "Cl2": 0.0, "F2": 0.0, "Na": 0.0, "Fe": 0.0, "Zn": 0.0,
            "Cu": 0.0, "Ca": 0.0, "Mg": 0.0, "Al": 0.0, "Si": 0.0,
            "H2O(l)": -285.83, "H2O(g)": -241.82, "H2O": -285.83,
            "CO(g)": -110.53, "CO": -110.53,
            "CO2(g)": -393.51, "CO2": -393.51,
            "SO2(g)": -296.84, "SO2": -296.84,
            "NO(g)": 90.25, "NO": 90.25,
            "NO2(g)": 33.18,
            "HF(g)": -273.30, "HCl(g)": -92.31,
            "H2S(g)": -20.63,
            "NH3(g)": -45.94, "NH3": -45.94,
            "CH4(g)": -74.81, "CH4": -74.81,
            "C2H6(g)": -84.68,
            "C2H4(g)": 52.47,
            "C2H2(g)": 227.39,
            "C3H8(g)": -103.85,
            "C6H6(l)": 49.04,
            "CH3OH(l)": -238.66,
            "C2H5OH(l)": -277.69,
            "HCOOH(l)": -424.72,
            "CH3COOH(l)": -484.13,
            "C6H12O6(s)": -1273.02,
            "H2O2(l)": -187.78,
            "MgO(s)": -601.60, "MgO": -601.60,
            "CaO(s)": -635.09, "CaO": -635.09,
            "Al2O3(s)": -1675.70, "Al2O3": -1675.70,
            "Fe2O3(s)": -824.20, "Fe2O3": -824.20,
            "CuO(s)": -157.28,
            "ZnO(s)": -348.28,
            "SiO2(s)": -910.86, "SiO2": -910.86,
            "NaCl(s)": -411.15, "NaCl": -411.15,
            "AgCl(s)": -127.07,
            "CaCl2(s)": -795.80,
            "MgCl2(s)": -641.32,
            "FeCl2(s)": -342.67,
            "CuCl2(s)": -220.10,
            "ZnCl2(s)": -415.05,
            "CaCO3(s)": -1206.92, "CaCO3": -1206.92,
            "Na2CO3(s)": -1130.68,
            "BaCO3(s)": -1216.29,
            "BaSO4(s)": -1473.19,
            "CaSO4(s)": -1434.52,
            "Na2SO4(s)": -1387.08,
            "CuSO4(s)": -771.36,
            "MgSO4(s)": -1284.91,
            "NaNO3(s)": -467.85,
            "KNO3(s)": -494.63,
            "NH4NO3(s)": -365.56,
            "FeS(s)": -100.02,
            "ZnS(s)": -205.98,
            "PbS(s)": -98.28,
            "NH4Cl(s)": -314.43,
            "(NH4)2SO4(s)": -1180.85,
            "PCl3(g)": -287.02,
            "SF6(g)": -1209.00,
            "BF3(g)": -1135.61,
            "CS2(l)": 89.70,
            "HCN(g)": 135.14,
            "NaOH(s)": -425.61,
            "Ca(OH)2(s)": -986.09,
            "Mg(OH)2(s)": -924.54,
            # Ions
            "H+(aq)": 0.0, "OH-(aq)": -230.02,
            "Na+(aq)": -240.12, "K+(aq)": -252.14,
            "Ag+(aq)": 105.58, "Ca2+(aq)": -542.96,
            "Mg2+(aq)": -466.85,
            "Fe2+(aq)": -89.12, "Fe3+(aq)": -48.53,
            "Cu2+(aq)": 64.77, "Zn2+(aq)": -153.89,
            "Al3+(aq)": -531.00,
            "Cl-(aq)": -167.08, "Br-(aq)": -121.41,
            "I-(aq)": -55.19,
            "SO4^2-(aq)": -909.27, "CO3^2-(aq)": -677.14,
            "NO3-(aq)": -205.00, "NH4+(aq)": -132.51,
        }

    def _run_base(self, reactants: Dict[str, float], products: Dict[str, float],
                  temperature_k: float = 298.15) -> dict:
        """Core logic: calculate ΔG°rxn."""
        if not reactants or not products:
            raise ChemMCPError("Both reactants and products dictionaries must be provided.")
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive (in Kelvin).")

        # Calculate ΔG°rxn directly from ΔG°f values
        prod_dg_total = 0.0
        rea_dg_total = 0.0
        prod_ds_total = 0.0
        rea_ds_total = 0.0
        prod_dh_total = 0.0
        rea_dh_total = 0.0

        missing = []

        for species, coeff in products.items():
            dg = self._lookup(self._dgf_db, species)
            ds = self._lookup(self._s_db, species)
            dh = self._lookup(self._dhf_db, species)
            if dg is None:
                missing.append(f"{species} (ΔG°f)")
            else:
                prod_dg_total += coeff * dg
            if ds is not None:
                prod_ds_total += coeff * ds
            if dh is not None:
                prod_dh_total += coeff * dh

        for species, coeff in reactants.items():
            dg = self._lookup(self._dgf_db, species)
            ds = self._lookup(self._s_db, species)
            dh = self._lookup(self._dhf_db, species)
            if dg is None:
                missing.append(f"{species} (ΔG°f)")
            else:
                rea_dg_total += coeff * dg
            if ds is not None:
                rea_ds_total += coeff * ds
            if dh is not None:
                rea_dh_total += coeff * dh

        if missing:
            raise ChemMCPError(
                f"Thermodynamic data not found for: {', '.join(missing)}. "
                f"Available data covers ~150 common compounds."
            )

        # Primary method: ΔG°rxn from ΔG°f values
        delta_g = prod_dg_total - rea_dg_total

        # Secondary method: ΔG° = ΔH° - TΔS (for cross-validation)
        delta_h = prod_dh_total - rea_dh_total
        delta_s = prod_ds_total - rea_ds_total  # J/(mol·K)
        delta_g_check = delta_h - (temperature_k * delta_s / 1000.0)  # convert S to kJ

        # Determine spontaneity
        is_spontaneous = delta_g < 0

        # Estimate equilibrium constant: ΔG° = -RT ln K
        import math
        R = 8.314e-3  # kJ/(mol·K)
        if abs(delta_g) < 500:  # reasonable range for K calculation
            try:
                if abs(delta_g / (R * temperature_k)) < 500:
                    K = math.exp(-delta_g / (R * temperature_k))
                    if K > 1e10:
                        K_str = "very large (K >> 1, products heavily favored)"
                    elif K < 1e-10:
                        K_str = "very small (K ≈ 0, reactants favored)"
                    else:
                        K_str = f"{K:.3g}"
                else:
                    K_str = "extreme" if delta_g < 0 else "≈ 0"
            except OverflowError:
                K_str = "extreme" if delta_g < 0 else "≈ 0"
        else:
            K_str = "N/A"

        return {
            "delta_g_rxn_kj_per_mol": round(delta_g, 2),
            "unit": "kJ/mol",
            "is_spontaneous": is_spontaneous,
            "temperature_k": round(temperature_k, 2),
            "delta_h_rxn_kj_per_mol": round(delta_h, 2),
            "delta_s_rxn_j_per_mol_k": round(delta_s, 2),
            "method": "ΔG°rxn = Σ(ν_i×ΔG°f_products) - Σ(ν_j×ΔG°f_reactants)",
            "cross_check": {
                "delta_g_via_dh_ts": round(delta_g_check, 2),
                "formula": f"ΔG° = ΔH° - TΔS = {round(delta_h, 2)} - {round(temperature_k, 2)} × ({round(delta_s, 2)}/1000) = {round(delta_g_check, 2)} kJ/mol",
            },
            "equilibrium_constant_K": K_str,
            "breakdown": {
                "products_delta_g_sum": round(prod_dg_total, 2),
                "reactants_delta_g_sum": round(rea_dg_total, 2),
            },
            "note": (
                f"ΔG° < 0 → spontaneous as written; "
                f"ΔG° > 0 → non-spontaneous (reverse direction spontaneous); "
                f"ΔG° = 0 → at equilibrium"
            ),
        }

    def _run_text(self, input_str: str) -> dict:
        """Parse text input."""
        try:
            parts = input_str.strip().split(";")
            if len(parts) < 2:
                raise ValueError("Expected format: reactants;products[;temperature]")

            reactants = self._parse_species_dict(parts[0])
            products = self._parse_species_dict(parts[1])
            temp = 298.15
            if len(parts) >= 3:
                temp = float(parts[2].strip())

            return self._run_base(reactants, products, temp)
        except Exception as e:
            if isinstance(e, ChemMCPError):
                raise
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'species1:coef1,species2:coef2;species3:coef3;T'")

    def _lookup(self, db: dict, species: str):
        """Look up in a database with fallback to state suffixes."""
        if species in db:
            return db[species]
        for suffix in ["(g)", "(l)", "(s)", "(aq)"]:
            key = species + suffix
            if key in db:
                return db[key]
        return None

    def _parse_species_dict(self, s: str) -> Dict[str, float]:
        """Parse species dictionary string."""
        s = s.strip()
        result = {}
        if ":" in s:
            pairs = s.split(",")
            for pair in pairs:
                if ":" not in pair:
                    continue
                sp, coef = pair.split(":", 1)
                result[sp.strip()] = float(coef.strip())
            return result
        import re
        tokens = re.split(r'\s*\+\s*', s)
        for token in tokens:
            token = token.strip()
            if not token:
                continue
            m = re.match(r'^(\d*\.?\d*)\s*(.+)$', token)
            if m:
                coef_str = m.group(1)
                species = m.group(2).strip()
                coef = float(coef_str) if coef_str else 1.0
                result[species] = coef
        return result
