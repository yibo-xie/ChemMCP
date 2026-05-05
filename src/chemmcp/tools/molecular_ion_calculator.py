import logging
import math
from typing import Optional, Dict, Any, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Common adduct masses (in Da, relative to neutral mass)
# Format: [adduct_key] = (mass_shift_Da, charge, description)
ADDUCTS = {
    # Positive mode
    "[M+H]+": (1.007276, 1, "Protonation — most common in ESI+"),
    "[M+Na]+": (22.989218, 1, "Sodium adduct — common for sugars, PEGs"),
    "[M+K]+": (38.963158, 1, "Potassium adduct"),
    "[M+NH4]+": (18.033823, 1, "Ammonium adduct — common for amines, lipids"),
    "[M+2H]2+": (0.503638, 2, "Doubly protonated — peptides, basic compounds"),
    "[M+H-H2O]+": (-17.002740, 1, "Dehydrated protonated — loss of water from [M+H]"),
    "[M+CH5O2]+": (51.997160, 1, "Methanol adduct + H — MeOH solvent effect"),
    "[M+ACN+H]+": (42.034398, 1, "Acetonitrile adduct + H — ACN solvent effect"),
    "[M+2Na]2+": (11.494609, 2, "Disodium adduct, 2+"),

    # Negative mode
    "[M-H]-": (-1.007276, -1, "Deprotonation — most common in ESI-"),
    "[M+Cl]-": (34.968853, -1, "Chloride adduct — common for acidic compounds"),
    "[M+FA-H]-": (44.997604, -1, "Formate adduct — formic acid modifier"),
    "[M+Ac-H]-": (59.013304, -1, "Acetate adduct — acetic acid modifier"),
    "[M-2H]2-": (-0.503638, -2, "Doubly deprotonated — acidic compounds"),
    "[M+CF3COO]-": (112.985589, -1, "Trifluoroacetate adduct — TFA contamination warning"),
    "[M+HCOO]-": (44.997204, -1, "Formate adduct (alternative)"),

    # Other
    "[M]+.": (0.000548, 1, "Molecular ion radical cation (EI)"),  # ~electron mass
}


# Common elemental isotopes for isotope pattern reference
ISOTOPE_DATA = {
    "C": {"exact_mass": 12.000000, "A1_abundance": 1.11, "A1_mass_diff": 1.003355},  # C-13
    "H": {"exact_mass": 1.007825, "A1_abundance": 0.016, "A1_mass_diff": 0.004019},
    "N": {"exact_mass": 14.003074, "A1_abundance": 0.37, "A1_mass_diff": 0.000000},  # N-15 same nominal mass
    "O": {"exact_mass": 15.994915, "A1_abundance": 0.04, "A1_mass_diff": 0.004034},
    "S": {"exact_mass": 31.972071, "A1_abundance": 4.42, "A1_mass_diff": 0.999388},
    "Cl": {"exact_mass": 34.968853, "A1_abundance": 32.50, "A1_mass_diff": 1.997050},
    "Br": {"exact_mass": 78.918336, "A1_abundance": 49.31, "A1_mass_diff": 1.997905},
}


@ChemMCPManager.register_tool
class MolecularIonCalculator(BaseTool):
    """
    分子离子峰质量计算工具（考虑加合物）。
    计算分子离子峰的精确质量、同位素分布、常见加合物m/z值，辅助质谱解析。
    """
    __version__ = "0.1.0"
    name = "MolecularIonCalculator"
    func_name = "calculate_molecular_ion"
    description = "Calculate exact molecular ion m/z values considering adducts, charge states, and isotope patterns for mass spectrometry analysis."
    implementation_description = (
        "Computes exact monoisotopic mass and m/z for a given molecular formula across multiple "
        "adduct types (ESI+, ESI-, EI, APCI). Includes comprehensive adduct database "
        "([M+H]⁺, [M+Na]⁺, [M+NH₄]⁺, [M-H]⁻, [M+Cl]⁻, etc.). "
        "Calculates theoretical isotope patterns using elemental compositions. "
        "Provides m/z tables for rapid MS peak identification."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Mass Spectrometry", "Molecular Ion", "Adduct", "Exact Mass", "Isotope Pattern", "MS", "m/z"]
    required_envs = []

    code_input_sig = [
        ("molecular_formula", "str", "N/A", "Molecular formula (e.g., 'C16H13NO2', 'C6H12O6')."),
        ("target_adducts", "list", "None", "List of adduct keys to calculate (e.g., ['[M+H]+', '[M+Na]+']). If None, calculate all common ones."),
        ("ionization_mode", "str", "'positive'", "Ionization mode: 'positive', 'negative', or 'both'."),
        ("charge_state", "int", "None", "Specific charge state override (e.g., 2 for 2+). None = auto-detect from adduct."),
        ("mass_tolerance_ppm", "float", "5.0", "Mass tolerance in ppm for matching/identification."),
        ("isotope_max_peaks", "int", "5", "Number of isotope peaks to calculate (M, M+1, M+2, ...)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Parameters like 'formula=C16H13NO2 mode=positive' or 'formula=C6H12O6 adducts=[M+H]+,[M+Na]+'."),
    ]

    output_sig = [
        ("ion_analysis", "dict", "Complete molecular ion analysis including exact mass, adduct m/z table, isotope pattern, and identification guidance."),
    ]

    examples = [
        {
            "code_input": {
                "molecular_formula": "C16H13NO2",
            },
            "text_input": {"input_params": "formula=C16H13NO2"},
            "output": {
                "ion_analysis": {"monoisotopic_mass": 251.094629, "adduct_mz_table": {...}}
            },
        },
        {
            "code_input": {
                "molecular_formula": "C6H12O6",
                "target_adducts": ["[M+H]+", "[M+Na]+", "[M-H]-"],
            },
            "text_input": {"input_params": "formula=C6H12O6 adducts=[M+H]+,[M+Na]+,[M-H]-"},
            "output": {"ion_analysis": {...}},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    # Elemental atomic weights (monoisotopic masses)
    ATOMIC_MASS = {
        "H": 1.00782503223, "He": 3.01602932009, "Li": 7.0160034366,
        "B": 10.811, "Be": 9.0121831865, "C": 12.0000000, "N": 14.00307400443,
        "O": 15.99491461957, "F": 18.99840316273, "Ne": 20.1797, "Na": 22.9897692820,
        "Mg": 23.985041697, "Al": 26.981538453, "Si": 27.976926535, "P": 30.973761998,
        "S": 31.9720711744, "Cl": 34.968852682, "Ar": 39.9623831225, "K": 38.96370648,
        "Ca": 39.962590863, "Sc": 44.95590828, "Ti": 47.94794198, "V": 50.94795125,
        "Cr": 51.94050623, "Mn": 54.93804391, "Fe": 55.9349375, "Co": 58.93319429,
        "Ni": 58.93334377, "Cu": 63.546, "Zn": 63.92914215, "Ga": 69.723, "Ge": 72.630,
        "As": 74.9215942, "Se": 78.971, "Br": 78.9183376, "Kr": 83.798, "Rb": 84.911789728,
        "Sr": 87.905635, "Y": 88.9058403, "Zr": 91.224, "Nb": 92.906373, "Mo": 95.95,
        "I": 126.9044719, "Xe": 131.293, "Cs": 132.90545196, "Ba": 137.327, "La": 138.90547,
    }

    def _parse_formula(self, formula: str) -> dict:
        """Parse molecular formula into element counts."""
        import re
        formula = formula.strip()
        elements = {}
        # Match element symbols (capital letter + optional lowercase) followed by optional number
        pattern = r"([A-Z][a-z]?)(\d*)"
        for match in re.finditer(pattern, formula):
            elem = match.group(1)
            count_str = match.group(2)
            count = int(count_str) if count_str else 1
            if elem:
                elements[elem] = elements.get(elem, 0) + count

        if not elements:
            raise ChemMCPError(f"Could not parse molecular formula '{formula}'.")
        return elements

    def _calc_monoisotopic_mass(self, elements: dict) -> float:
        """Calculate exact monoisotopic mass."""
        mass = 0.0
        unknown = []
        for elem, count in elements.items():
            if elem in self.ATOMIC_MASS:
                mass += self.ATOMIC_MASS[elem] * count
            else:
                unknown.append(elem)

        if unknown:
            raise ChemMCPError(f"Unknown element(s): {', '.join(unknown)}. Supported: {', '.join(sorted(set(self.ATOMIC_MASS.keys()) & set(elements.keys())))}")

        return mass

    def _calc_isotope_pattern(self, elements: dict, n_peaks: int = 5) -> list:
        """
        Calculate approximate isotope distribution.
        Simplified model: M+1 from C-13, M+2 from combinations of C-13/C-13, O-18, S-34, etc.
        """
        # For simplicity, use polynomial expansion approximation
        # Focus on major contributors to M+1 and M+2

        n_C = elements.get("C", 0)
        n_H = elements.get("H", 0)
        n_N = elements.get("N", 0)
        n_O = elements.get("O", 0)
        n_S = elements.get("S", 0)
        n_Cl = elements.get("Cl", 0)
        n_Br = elements.get("Br", 0)

        # M+1 probability (mainly C-13)
        p_M1 = (n_C * 0.0111 +
                n_N * 0.0037 +
                n_O * 0.0004 +
                n_S * 0.0080 +
                n_H * 0.00015 +
                n_Cl * 0.0 +
                n_Br * 0.0)

        # M+2 probability (C-13_2, O-18, S-34, Cl-37, Br-81)
        p_M2_from_C13 = (n_C * (n_C - 1) / 2) * (0.0111 ** 2)
        p_M2_from_O18 = n_O * 0.0020
        p_M2_from_S34 = n_S * 0.0445
        p_M2_from_Cl37 = n_Cl * 0.3250 if n_Cl > 0 else 0
        p_M2_from_Br81 = n_Br * 0.4969 if n_Br > 0 else 0
        p_M2 = p_M2_from_C13 + p_M2_from_O18 + p_M2_from_S34 + p_M2_from_Cl37 + p_M2_from_Br81

        # M+3 (mainly from halogens or large molecules)
        p_M3 = 0
        if n_Cl >= 2:
            p_M3 += (n_Cl * (n_Cl - 1) / 2) * (0.325 ** 2)
        if n_Br >= 2:
            p_M3 += (n_Br * (n_Br - 1) / 2) * (0.4969 ** 2)
        if n_Cl >= 1 and n_Br >= 1:
            p_M3 += n_Cl * n_Br * 0.325 * 0.4969

        # M+4 (for di-halogenated compounds)
        p_M4 = 0
        if n_Cl >= 2:
            # Skip complex calculation; just note it exists
            pass
        if n_Br >= 2:
            p_M4 += (n_Br * (n_Br - 1) * (n_Br - 2) / 6) * (0.4969 ** 3)

        base_mass = self._calc_monoisotopic_mass(elements)
        pattern = [
            {"peak": "M", "mz": round(base_mass, 6), "relative_abundance_pct": 100.0, "mass_shift_Da": 0},
        ]
        if n_peaks > 1 and p_M1 > 0.001:
            pattern.append({"peak": "M+1", "mz": round(base_mass + 1.003355, 6),
                           "relative_abundance_pct": round(p_M1 * 100, 2), "mass_shift_Da": 1.003355})
        if n_peaks > 2 and p_M2 > 0.001:
            pattern.append({"peak": "M+2", "mz": round(base_mass + 2.004671, 6),
                           "relative_abundance_pct": round(p_M2 * 100, 2), "mass_shift_Da": 2.004671})
        if n_peaks > 3 and p_M3 > 0.01:
            pattern.append({"peak": "M+3", "mz": round(base_mass + 3.006019, 6),
                           "relative_abundance_pct": round(p_M3 * 100, 2), "mass_shift_Da": 3.006019})
        if n_peaks > 4 and p_M4 > 0.1:
            pattern.append({"peak": "M+4", "mz": round(base_mass + 4.007820, 6),
                           "relative_abundance_pct": round(p_M4 * 100, 2), "mass_shift_Da": 4.007820})

        return pattern

    def _calc_adduct_mz(self, mono_mass: float, target_adducts: Optional[List[str]],
                         mode: str) -> list:
        """Calculate m/z for each relevant adduct."""
        results = []
        adduct_keys = list(ADDUCTS.keys())

        if target_adducts:
            adduct_keys = [a for a in adduct_keys if a in target_adducts]

        if mode == "positive":
            adduct_keys = [a for a in adduct_keys if a.endswith("+") or "+" in a]
        elif mode == "negative":
            adduct_keys = [a for a in adduct_keys if a.endswith("-") or "-" in a]

        for key in sorted(adduct_keys):
            shift, charge, desc = ADDUCTS[key]
            mz = (mono_mass + shift) / abs(charge)
            results.append({
                "adduct": key,
                "charge": charge,
                "mz_exact": round(mz, 4),
                "mass_shift_Da": shift,
                "description": desc,
            })

        return results

    def _identify_elemental_composition(self, elements: dict) -> dict:
        """Summarize composition characteristics."""
        total_atoms = sum(elements.values())
        heteroatoms = {k: v for k, v in elements.items() if k not in ("C", "H")}
        halogens = {k: v for k, v in elements.items() if k in ("Cl", "Br", "F", "I")}

        return {
            "total_atoms": total_atoms,
            "carbon_count": elements.get("C", 0),
            "hydrogen_count": elements.get("H", 0),
            "heteroatom_types": list(heteroatoms.keys()),
            "halogen_present": bool(halogens),
            "halogen_details": dict(halogens) if halogens else None,
            "degree_of_unsaturation": self._calc_doubling_equivalent(elements),
            "nitrogen_rule_check": self._check_nitrogen_rule(mono_mass=self._calc_monoisotopic_mass(elements), elements=elements),
        }

    @staticmethod
    def _calc_doubling_equivalent(elements: dict) -> float:
        """Calculate degree of unsaturation (double bond equivalents)."""
        C = elements.get("C", 0)
        H = elements.get("H", 0)
        N = elements.get("N", 0)
        X = elements.get("F", 0) + elements.get("Cl", 0) + elements.get("Br", 0) + elements.get("I", 0)
        P = elements.get("P", 0)
        return C + 1 - (H + X - N)/2.0 + (N + P)/2.0

    @staticmethod
    def _check_nitrogen_rule(mono_mass: float, elements: dict) -> dict:
        """Nitrogen rule: odd N → odd MW; even N → even MW."""
        n_N = elements.get("N", 0)
        mw_is_even = int(round(mono_mass)) % 2 == 0
        n_even = n_N % 2 == 0
        consistent = (mw_is_even == n_even)
        return {
            "nitrogen_count": n_N,
            "nominal_mass_parity": "even" if mw_is_even else "odd",
            "rule_satisfied": consistent,
            "note": "✅ Nitrogen rule satisfied" if consistent else "⚠️ Check formula — nitrogen rule violation",
        }

    def _run_base(self, molecular_formula: str,
                  target_adducts: Optional[List[str]] = None,
                  ionization_mode: str = "positive",
                  charge_state: Optional[int] = None,
                  mass_tolerance_ppm: float = 5.0,
                  isotope_max_peaks: int = 5) -> dict:
        """Core logic."""

        # Parse formula
        elements = self._parse_formula(molecular_formula)

        # Monoisotopic mass
        mono_mass = self._calc_monoisotopic_mass(elements)

        # Adduct m/z table
        adduct_results = self._calc_adduct_mz(mono_mass, target_adducts, ionization_mode)

        # Isotope pattern
        isotope_pattern = self._calc_isotope_pattern(elements, isotope_max_peaks)

        # Composition summary
        composition = self._identify_elemental_composition(elements)

        result = {
            "ion_analysis": {
                "input": {
                    "molecular_formula": molecular_formula,
                    "ionization_mode": ionization_mode,
                    "mass_tolerance_ppm": mass_tolerance_ppm,
                },
                "monoisotopic_mass": {
                    "exact_mass_Da": round(mono_mass, 6),
                    "nominal_mass_Da": int(round(mono_mass)),
                },
                "adduct_mz_table": adduct_results,
                "isotope_pattern": isotope_pattern,
                "elemental_composition": composition,
                "ms_identification_tips": self._get_ms_tips(
                    molecular_formula, elements, mono_mass, adduct_results, composition),
            }
        }
        return result

    def _get_ms_tips(self, formula: str, elements: dict, mono_mass: float,
                      adducts: list, composition: dict) -> List[str]:
        tips = []

        # Primary ion to look for
        pos_adducts = [a for a in adducts if a["charge"] > 0]
        neg_adducts = [a for a in adducts if a["charge"] < 0]

        if pos_adducts:
            best_pos = min(pos_adducts, key=lambda x: x["mz_exact"])
            tips.append(f"🔍 ESI+: Look for m/z {best_pos['mz_exact']:.4f} ({best_pos['adduct']}) as primary signal")

        if neg_adducts:
            best_neg = min(neg_adducts, key=lambda x: abs(x["mz_exact"]))
            tips.append(f"🔍 ESI-: Look for m/z {best_neg['mz_exact']:.4f} ({best_neg['adduct']}) as primary signal")

        # Halogen pattern tip
        halogens = composition.get("halogen_details")
        if halogens:
            tips.append(f"⚛️ Halogen(s) present: {halogens} — check characteristic isotope pattern")
            if halogens.get("Cl", 0) == 1:
                tips.append("   ¹Cl: expect ~3:1 M:M+2 ratio")
            if halogens.get("Cl", 0) >= 2:
                tips.append(f"   ²Cl{halogens['Cl']}: expect characteristic 9:6:1 (or similar) M:M+1:M+2 pattern")
            if halogens.get("Br", 0) == 1:
                tips.append("   ⁸¹Br: expect ~1:1 M:M+2 ratio")

        # DBE tip
        dbe = composition.get("degree_of_unsaturation", 0)
        if dbe >= 4:
            tips.append(f"📐 High DBE ({dbe:.0f}) suggests aromatic rings or multiple π-systems")

        # Sodium adduct tip for certain compound types
        n_O = elements.get("O", 0)
        if n_O >= 4:
            tips.append("💡 Oxygen-rich compound — also check for [M+Na]⁺ and [M+NH₄]⁺ adducts")

        return tips[:7]

    def _run_text(self, input_params: str) -> dict:
        kwargs = {}
        parts = input_params.strip().split()
        for part in parts:
            if "=" in part:
                k, v = part.split("=", 1)
                if k == "formula":
                    kwargs["molecular_formula"] = v
                elif k == "mode":
                    kwargs["ionization_mode"] = v
                elif k == "adducts":
                    kwargs["target_adducts"] = [a.strip() for a in v.split(",")]
                elif k == "charge":
                    kwargs["charge_state"] = int(v)
                elif k == "ppm":
                    kwargs["mass_tolerance_ppm"] = float(v)
                elif k == "iso_peaks":
                    kwargs["isotope_max_peaks"] = int(v)
        return self._run_base(**kwargs)
