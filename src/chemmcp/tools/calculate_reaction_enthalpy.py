import logging
from typing import Dict, List, Optional, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CalculateReactionEnthalpy(BaseTool):
    """
    用 Hess 定律计算反应焓变 ΔH°rxn。
    ΔH°rxn = Σ(ν_i × ΔH°f(产物)) - Σ(ν_j × ΔH°f(反应物))
    """
    __version__ = "0.1.0"
    name = "CalculateReactionEnthalpy"
    func_name = "calculate_reaction_enthalpy"
    description = "Calculate reaction enthalpy change (ΔH°rxn) using Hess's Law from standard enthalpies of formation."
    implementation_description = "Applies Hess's Law: ΔH°rxn = Σ(n_i × ΔH°f_products) - Σ(m_j × ΔH°f_reactants). Uses built-in thermodynamic database of standard enthalpies of formation."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Hess Law", "Enthalpy", "Thermochemistry", "Reaction Heat"]
    required_envs = []

    code_input_sig = [
        ("reactants", "dict", "N/A", "Reactants as {species: stoichiometric_coefficient}, e.g., {'C3H8': 5, 'O2': 2}. Coefficients can be int or float."),
        ("products", "dict", "N/A", "Products as {species: stoichiometric_coefficient}, e.g., {'CO2': 3, 'H2O': 4}."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Semicolon-separated: reactants;products. E.g., 'C3H8:5,O2:2;CO2:3,H2O:4' or 'CH4+2O2;CO2+2H2O'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with delta_h_rxn (kJ/mol), breakdown by species, hess_law_application, and reaction_type (exothermic/endothermic)."),
    ]

    examples = [
        {
            "code_input": {
                "reactants": {"C3H8": 1, "O2": 5},
                "products": {"CO2": 3, "H2O": 4},
            },
            "text_input": {
                "input_str": "C3H8:1,O2:5;CO2:3,H2O:4",
            },
            "output": {
                "result": {
                    "delta_h_rxn_kj_per_mol": -2043.99,
                    "unit": "kJ/mol",
                    "reaction_type": "exothermic",
                    "hess_law_application": "ΔH°rxn = [3×(-393.51) + 4×(-285.83)] - [1×(-103.85) + 5×0] = -2043.99 kJ/mol",
                    "breakdown": {
                        "products_total": -2433.35,
                        "reactants_total": -103.85,
                        "details": {"CO2": (-393.51, 3), "H2O": (-285.83, 4), "C3H8": (-103.85, 1), "O2": (0.0, 5)},
                    },
                }
            },
        },
        {
            "code_input": {
                "reactants": {"N2": 1, "H2": 3},
                "products": {"NH3": 2},
            },
            "text_input": {
                "input_str": "N2:1,H2:3;NH3:2",
            },
            "output": {
                "result": {
                    "delta_h_rxn_kj_per_mol": -91.88,
                    "unit": "kJ/mol",
                    "reaction_type": "exothermic",
                    "hess_law_application": "ΔH°rxn = [2×(-45.94)] - [1×0 + 3×0] = -91.88 kJ/mol (Haber process)",
                    "breakdown": {
                        "products_total": -91.88,
                        "reactants_total": 0.0,
                        "details": {"NH3": (-45.94, 2), "N2": (0.0, 1), "H2": (0.0, 3)},
                    },
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize standard enthalpy of formation database for Hess's Law calculations."""
        # Same database as GetStandardEnthalpy — key thermodynamic data
        # ΔH°f in kJ/mol at 298 K
        self._dhf_db = {
            # Elements (0)
            "H2": 0.0, "O2": 0.0, "N2": 0.0, "C": 0.0, "S": 0.0,
            "Cl2": 0.0, "Br2": 0.0, "I2": 0.0, "F2": 0.0,
            "Na": 0.0, "Fe": 0.0, "Zn": 0.0, "Cu": 0.0, "Ag": 0.0,
            "Ca": 0.0, "Mg": 0.0, "Al": 0.0, "Si": 0.0,
            # Common compounds
            "H2O(l)": -285.83, "H2O(g)": -241.82, "H2O": -285.83,
            "CO(g)": -110.53, "CO": -110.53,
            "CO2(g)": -393.51, "CO2": -393.51,
            "SO2(g)": -296.84, "SO2": -296.84,
            "SO3(g)": -395.72,
            "NO(g)": 90.25, "NO": 90.25,
            "NO2(g)": 33.18,
            "N2O(g)": 82.05,
            "HF(g)": -273.30,
            "HCl(g)": -92.31,
            "HBr(g)": -36.40,
            "HI(g)": 26.48,
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
            "HCHO(g)": -108.57,
            "HCOOH(l)": -424.72,
            "CH3COOH(l)": -484.13,
            "C6H12O6(s)": -1273.02,
            "C12H22O11(s)": -2226.09,
            "H2O2(l)": -187.78,
            "H2SO4(l)": -814.00,
            "HNO3(l)": -174.10,
            "H3PO4(s)": -1279.01,
            # Oxides
            "MgO(s)": -601.60, "MgO": -601.60,
            "CaO(s)": -635.09, "CaO": -635.09,
            "Al2O3(s)": -1675.70, "Al2O3": -1675.70,
            "Fe2O3(s)": -824.20, "Fe2O3": -824.20,
            "Fe3O4(s)": -1118.38,
            "FeO(s)": -272.00,
            "CuO(s)": -157.28,
            "Cu2O(s)": -168.60,
            "ZnO(s)": -348.28,
            "SiO2(s)": -910.86, "SiO2": -910.86,
            "Cr2O3(s)": -1139.70,
            "MnO2(s)": -520.03,
            "PbO(s)": -217.32,
            "PbO2(s)": -277.40,
            "TiO2(s)": -944.00,
            "BaO(s)": -548.10,
            # Acids / Bases
            "NaOH(s)": -425.61,
            "KOH(s)": -424.72,
            "Ca(OH)2(s)": -986.09,
            "Al(OH)3(s)": -1277.00,
            "Fe(OH)3(s)": -823.00,
            "Mg(OH)2(s)": -924.54,
            # Halides
            "NaCl(s)": -411.15, "NaCl": -411.15,
            "KCl(s)": -436.75,
            "AgCl(s)": -127.07,
            "AgBr(s)": -100.37,
            "AgI(s)": -61.84,
            "CaCl2(s)": -795.80,
            "MgCl2(s)": -641.32,
            "AlCl3(s)": -704.20,
            "FeCl2(s)": -342.67,
            "FeCl3(s)": -399.49,
            "CuCl2(s)": -220.10,
            "ZnCl2(s)": -415.05,
            "PbCl2(s)": -359.41,
            # Carbonates
            "CaCO3(s)": -1206.92, "CaCO3": -1206.92,
            "Na2CO3(s)": -1130.68,
            "NaHCO3(s)": -950.81,
            "BaCO3(s)": -1216.29,
            "MgCO3(s)": -1095.79,
            # Sulfates
            "CaSO4(s)": -1434.52,
            "BaSO4(s)": -1473.19,
            "Na2SO4(s)": -1387.08,
            "CuSO4(s)": -771.36,
            "FeSO4(s)": -828.40,
            "MgSO4(s)": -1284.91,
            "ZnSO4(s)": -982.82,
            # Nitrates
            "NaNO3(s)": -467.85,
            "KNO3(s)": -494.63,
            "Ca(NO3)2(s)": -938.21,
            "NH4NO3(s)": -365.56,
            # Sulfides
            "FeS(s)": -100.02,
            "ZnS(s)": -205.98,
            "PbS(s)": -98.28,
            "CuS(s)": -53.06,
            "CdS(s)": -161.90,
            # Ammonium salts
            "NH4Cl(s)": -314.43,
            "(NH4)2SO4(s)": -1180.85,
            # Other
            "PCl3(g)": -287.02,
            "PCl5(g)": -374.89,
            "SF6(g)": -1209.00,
            "BF3(g)": -1135.61,
            "CS2(l)": 89.70,
            "HCN(g)": 135.14,
            "PH3(g)": 5.44,
            "Na2O(s)": -414.22,
            "Na2O2(s)": -510.87,
            "P4O10(s)": -2984.03,
            "O3(g)": 142.67,
            # Ions (aqueous)
            "H+(aq)": 0.0, "OH-(aq)": -230.02,
            "Na+(aq)": -240.12, "K+(aq)": -252.14,
            "Ag+(aq)": 105.58, "Ca2+(aq)": -542.96,
            "Mg2+(aq)": -466.85, "Ba2+(aq)": -537.64,
            "Fe2+(aq)": -89.12, "Fe3+(aq)": -48.53,
            "Cu2+(aq)": 64.77, "Zn2+(aq)": -153.89,
            "Al3+(aq)": -531.00,
            "Cl-(aq)": -167.08, "Br-(aq)": -121.41,
            "I-(aq)": -55.19,
            "SO4^2-(aq)": -909.27, "CO3^2-(aq)": -677.14,
            "NO3-(aq)": -205.00, "NH4+(aq)": -132.51,
        }

    def _run_base(self, reactants: Dict[str, float], products: Dict[str, float]) -> dict:
        """Core logic: calculate ΔH°rxn using Hess's Law."""
        if not reactants or not products:
            raise ChemMCPError("Both reactants and products dictionaries must be provided.")

        # Calculate sum of ΔH°f for products
        prod_total = 0.0
        prod_details = {}
        missing_prod = []

        for species, coeff in products.items():
            dhf = self._lookup_dhf(species)
            if dhf is None:
                missing_prod.append(species)
            else:
                contribution = coeff * dhf
                prod_total += contribution
                prod_details[species] = {"dhf_kj_mol": dhf, "coeff": coeff, "contribution": round(contribution, 2)}

        # Calculate sum of ΔH°f for reactants
        rea_total = 0.0
        rea_details = {}
        missing_rea = []

        for species, coeff in reactants.items():
            dhf = self._lookup_dhf(species)
            if dhf is None:
                missing_rea.append(species)
            else:
                contribution = coeff * dhf
                rea_total += contribution
                rea_details[species] = {"dhf_kj_mol": dhf, "coeff": coeff, "contribution": round(contribution, 2)}

        # Check for missing data
        all_missing = missing_prod + missing_rea
        if all_missing:
            raise ChemMCPError(
                f"Standard enthalpy of formation not found for: {', '.join(all_missing)}. "
                f"Please check species names against available database entries."
            )

        # Hess's Law: ΔH°rxn = Σ(ΔH°f_products) - Σ(ΔH°f_reactants)
        delta_h = prod_total - rea_total

        # Determine reaction type
        if delta_h < 0:
            rxn_type = "exothermic (releases heat)"
        elif delta_h > 0:
            rxn_type = "endothermic (absorbs heat)"
        else:
            rxn_type = "thermoneutral (no heat change)"

        # Build detailed equation string
        parts = [f"{coeff}×({dhf})" if coeff != 1 else f"({dhf})"
                  for species, (dhf, coeff) in {**{s: (d["dhf_kj_mol"], d["coeff"]) for s, d in prod_details.items()},
                                                       **{s: (d["dhf_kj_mol"], d["coeff"]) for s, d in rea_details.items()}}.items()]
        # Simpler readable format
        prod_str = " + ".join([f"{int(c) if c == int(c) else c}×ΔH°f({s})" if c != 1 else f"ΔH°f({s})" for s, c in products.items()])
        rea_str = " + ".join([f"{int(c) if c == int(c) else c}×ΔH°f({s})" if c != 1 else f"ΔH°f({s})" for s, c in reactants.items()])
        hess_eq = f"ΔH°rxn = [{prod_str}] - [{rea_str}] = {round(prod_total, 2)} - ({round(rea_total, 2)}) = {round(delta_h, 2)} kJ/mol"

        return {
            "delta_h_rxn_kj_per_mol": round(delta_h, 2),
            "unit": "kJ/mol",
            "reaction_type": rxn_type,
            "hess_law_application": hess_eq,
            "breakdown": {
                "products_total": round(prod_total, 2),
                "reactants_total": round(rea_total, 2),
                "product_details": prod_details,
                "reactant_details": rea_details,
            },
            "temperature_k": 298.15,
            "note": "Calculated at standard conditions (298.15 K, 1 bar). Negative value → exothermic; Positive → endothermic.",
        }

    def _run_text(self, input_str: str) -> dict:
        """Parse text input: 'reactants;products' format."""
        try:
            parts = input_str.strip().split(";")
            if len(parts) < 2:
                raise ValueError("Expected 'reactants;products' format")

            reactants = self._parse_species_dict(parts[0])
            products = self._parse_species_dict(parts[1])
            return self._run_base(reactants, products)
        except Exception as e:
            if isinstance(e, ChemMCPError):
                raise
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Expected format: 'species1:coef1,species2:coef2;species3:coef3,species4:coef4'")

    def _parse_species_dict(self, s: str) -> Dict[str, float]:
        """Parse 'A:1,B:2' or 'A+2B' style into dict."""
        s = s.strip()
        result = {}

        # Try colon-separated first: A:1,B:2
        if ":" in s:
            pairs = s.split(",")
            for pair in pairs:
                pair = pair.strip()
                if ":" not in pair:
                    continue
                sp, coef = pair.split(":", 1)
                result[sp.strip()] = float(coef.strip())
            return result

        # Try plus-separated: A + 2B + 3C
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

    def _lookup_dhf(self, species: str) -> Optional[float]:
        """Look up ΔH°f for a species."""
        # Direct match
        if species in self._dhf_db:
            return self._dhf_db[species]
        # Try with common state suffixes
        for suffix in ["(g)", "(l)", "(s)", "(aq)"]:
            key = species + suffix
            if key in self._dhf_db:
                return self._dhf_db[key]
        return None
