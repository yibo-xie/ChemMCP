import logging
import re
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MolecularIonCalculator(BaseTool):
    """
    分子离子峰计算工具。
    计算分子离子峰 m/z、同位素峰（M+1, M+2）强度比及同位素分布模式。
    """
    __version__ = "0.1.0"
    name = "MolecularIonCalculator"
    func_name = "calculate_molecular_ion"
    description = "Calculate the molecular ion peak (m/z), exact mass, isotope pattern (M, M+1, M+2 intensities), and isotope distribution for a given molecular formula or SMILES."
    implementation_description = "Calculates exact monoisotopic mass from atomic weights, predicts M+1 intensity from ¹³C/²H/¹⁵N/¹⁷O/³³S contributions, M+2 from ¹⁸O/³⁴S/³⁷Cl/⁸¹Br. Includes characteristic isotope pattern recognition for Cl, Br, S, Si."
    oss_dependencies = [
        ("RDKit", "https://www.rdkit.org/", "BSD-3-Clause"),
        ("IUPAC atomic weights", "IUPAC Technical Report 2016", None),
    ]
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Mass Spectrometry", "Molecular Ion", "Isotope Pattern", "Exact Mass", "m/z", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("formula", "str", "N/A", "Molecular formula (e.g., 'C6H12O6') or SMILES string."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Molecular formula or SMILES. Example: 'C6H12O6' or 'CCO'"),
    ]

    output_sig = [
        ("ion_data", "dict", "Complete molecular ion data including exact mass, m/z, isotope pattern, and pattern recognition."),
    ]

    examples = [
        {
            "code_input": {"formula": "C6H12O6"},
            "text_input": {"input_params": "C6H12O6"},
            "output": {
                "ion_data": {
                    "formula": "C₆H₁₂O₆",
                    "exact_mass": 180.0634,
                    "mz_M": 180,
                    "M1_ratio": 6.7,
                    "M2_ratio": 0.24,
                }
            },
        },
        {
            "code_input": {"formula": "CCl4"},
            "text_input": {"input_params": "CCl4"},
            "output": {
                "ion_data": {
                    "formula": "CCl₄",
                    "mz_M": 152,
                    "pattern_type": "chlorine isotope cluster (multiple Cl)",
                }
            },
        },
    ]

    # ========== ISOTOPIC DATA ==========
    # Format: {element: [(mass, abundance), ...] sorted by abundance descending}
    _ISOTOPES = {
        "H":  [(1.007825, 100.0), (2.014102, 0.015)],         # H, D
        "C":  [(12.000000, 100.0), (13.003355, 1.07)],          # ¹²C, ¹³C
        "N":  [(14.003074, 100.0), (15.000109, 0.37)],          # ¹⁴N, ¹⁵N
        "O":  [(15.994915, 100.0), (16.999132, 0.038), (17.999160, 0.205)],  # ¹⁶O, ¹⁷O, ¹⁸O
        "F":  [(18.998403, 100.0)],                               # ¹⁹F only
        "Na": [(22.989770, 100.0)],                               # ²³Na only
        "Si": [(27.976927, 100.0), (28.976495, 5.08), (29.973770, 3.35)],  # ²⁸Si, ²⁹Si, ³⁰Si
        "P":  [(30.973762, 100.0)],                               # ³¹P only
        "S":  [(31.972071, 100.0), (32.971458, 0.75), (33.967867, 4.22)],  # ³²S, ³³S, ³⁴S
        "Cl": [(34.968853, 100.0), (36.965903, 32.5)],           # ³⁵Cl, ³⁷Cl (~3:1)
        "Br": [(78.918337, 100.0), (80.916291, 97.28)],          # ⁷⁹Br, ⁸¹Br (~1:1)
        "I":  [(126.904473, 100.0)],                              # ¹²⁷I only
        "B":  [(11.009306, 100.0), (10.012937, 19.9)],           # ¹¹B, ¹⁰B
        "Al": [(26.981538, 100.0)],                               # ²⁷Al only
        "Ca": [(39.962591, 100.0), (41.958618, 0.647)],          # ⁴⁰Ca, ⁴²Ca (simplified)
        "Fe": [(55.934939, 100.0), (56.935396, 2.119), (57.933276, 2.245)],  # ⁵⁶Fe, ⁵⁷Fe, ⁵⁸Fe
    }

    # Per-atom M+1 contribution (% relative to M peak)
    _M1_CONTRIBUTIONS = {
        "C": 1.07,     # ¹³C: ~1.1% per C atom
        "H": 0.015,    # ²D: negligible but included
        "N": 0.37,     # ¹⁵N: ~0.37% per N atom
        "O": 0.038,    # ¹⁷O: small contribution to M+1
        "S": 0.75,     # ³³S contributes to M+1
        "Si": 5.08,    # ²⁹Si significant!
        "Fe": 2.119,   # ⁵⁷Fe
    }

    # Per-atom M+2 contribution (% relative to M peak)
    _M2_CONTRIBUTIONS = {
        "O": 0.205,    # ¹⁸O: main M+2 contributor from oxygen
        "S": 4.22,     # ³⁴S: significant M+2 from sulfur
        "Si": 3.35,    # ³⁰Si: M+2 from silicon
        "Cl": 32.5,    # ³⁷Cl: HUGE M+2 (dominates!)
        "Br": 97.28,   # ⁸¹Br: ~equal M and M+2!
        "Fe": 2.245,   # ⁵⁸Fe
    }

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize RDKit if available."""
        self._rdkit_available = False
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors
            self.Chem = Chem
            self.Descriptors = Descriptors
            self._rdkit_available = True
        except ImportError:
            logger.warning("RDKit not available for MolecularIonCalculator")

    def _parse_formula(self, formula_or_smiles: str) -> dict:
        """Parse molecular formula into element counts."""
        s = formula_or_smiles.strip()

        # Try RDKit first if it looks like SMILES
        if self._rdkit_available and ('=' in s or '(' in s or '[' in s or len(s) > 10 or s[0].isupper() and s[0].islower()):
            try:
                mol = self.Chem.MolFromSmiles(s)
                if mol:
                    rd_formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
                    return self._parse_molecular_formula(rd_formula)
            except Exception:
                pass

        return self._parse_molecular_formula(s)

    def _parse_molecular_formula(self, formula: str) -> dict:
        """Parse molecular formula like C6H12O6 into element counts."""
        elements = {}
        # Handle nested groups like (CH3)2 - simplified parsing
        expanded = formula

        # Expand parenthesized groups (basic support)
        while '(' in expanded:
            match = re.search(r'\(([A-Z0-9a-z]+)\)(\d*)', expanded)
            if not match:
                break
            group_content = match.group(1)
            multiplier = int(match.group(2)) if match.group(2) else 1
            # Parse elements in group and multiply
            sub_elements = {}
            elem_matches = re.findall(r'([A-Z][a-z]?)(\d*)', group_content)
            for elem, cnt in elem_matches:
                if elem:
                    sub_elements[elem] = sub_elements.get(elem, 0) + int(cnt) * multiplier
            # Add to elements dict
            for elem, cnt in sub_elements.items():
                elements[elem] = elements.get(elem, 0) + cnt
            expanded = expanded[:match.start()] + expanded[match.end():]

        # Parse remaining elements
        matches = re.findall(r'([A-Z][a-z]?)(\d*)', expanded)
        for elem, count_str in matches:
            if not elem:
                continue
            count = int(count_str) if count_str else 1
            elements[elem] = elements.get(elem, 0) + count

        return {"formula": formula, "elements": elements}

    def _calculate_exact_mass(self, elements: dict) -> float:
        """Calculate monoisotopic exact mass using lightest isotopes."""
        mass = 0.0
        for elem, count in elements.items():
            if elem in self._ISOTOPES:
                mass += self._ISOTOPES[elem][0][0] * count  # [0][0] = lightest isotope mass
            else:
                logger.warning(f"Unknown element '{elem}', using approximate mass")
                # Approximate mass fallback
                approx_masses = {"H": 1.00783, "C": 12.0, "N": 14.003, "O": 16.0, "P": 31.0}
                mass += approx_masses.get(elem, 1.0) * count
        return round(mass, 4)

    def _calculate_isotope_pattern(self, elements: dict) -> dict:
        """Calculate M+1 and M+2 relative intensities."""
        m1_total = 0.0
        m2_total = 0.0

        for elem, count in elements.items():
            # M+1 contributions
            if elem in self._M1_CONTRIBUTIONS:
                m1_total += self._M1_CONTRIBUTIONS[elem] * count

            # M+2 contributions
            if elem in self._M2_CONTRIBUTIONS:
                m2_total += self._M2_CONTRIBUTIONS[elem] * count

        # Also account for combinations of two M+1 contributors giving M+2
        # e.g., two ¹³C atoms can contribute to M+2
        c_count = elements.get("C", 0)
        if c_count >= 2:
            # Probability of two 13C atoms: C(n,2) * (0.0107)^2 * 100
            import math
            m2_from_2C = math.comb(c_count, 2) * (0.0107 ** 2) * 100
            m2_total += m2_from_2C

        return {
            "M_plus_1_percent": round(m1_total, 2),
            "M_plus_2_percent": round(m2_total, 2),
            "M_plus_1_relative_to_M": f"{round(m1_total / 100, 4)} : 1" if m1_total > 0 else "negligible",
            "M_plus_2_relative_to_M": f"{round(m2_total / 100, 4)} : 1" if m2_total > 0 else "negligible",
        }

    def _recognize_pattern(self, elements: dict) -> str:
        """Recognize characteristic isotope patterns."""
        cl_count = elements.get("Cl", 0)
        br_count = elements.get("Br", 0)
        s_count = elements.get("S", 0)
        si_count = elements.get("Si", 0)

        patterns = []

        if cl_count == 1:
            patterns.append("Single Cl: M:M+2 ≈ 3:1 ratio (characteristic doublet)")
        elif cl_count == 2:
            patterns.append("Two Cl: M:M+2:M+4 ≈ 9:6:1 ratio (triplet)")
        elif cl_count == 3:
            patterns.append("Three Cl: M:M+2:M+4:M+6 ≈ 27:27:9:1 ratio")

        if br_count == 1:
            patterns.append("Single Br: M:M+2 ≈ 1:1 ratio (characteristic doublet of nearly equal height)")
        elif br_count == 2:
            patterns.append("Two Br: M:M+2:M+4 ≈ 1:2:1 ratio (triplet)")

        if cl_count > 0 and br_count > 0:
            patterns.append(f"Mixed Cl({cl_count})+Br({br_count}): complex combined isotope pattern")

        if s_count == 1:
            patterns.append("Single S: small M+2 peak (~4.2%) from ³⁴S")
        elif s_count >= 2:
            patterns.append(f"{s_count} S atoms: noticeable M+2 from ³⁴S combinations")

        if si_count > 0:
            patterns.append(f"Si present: enhanced M+1 ({si_count * 5.08}%) and M+2 ({si_count * 3.35}%) from Si isotopes")

        if cl_count == 0 and br_count == 0 and s_count == 0:
            patterns.append("No characteristic heavy isotope pattern (C/H/N/O/F/P compound): simple M+1 from ¹³C dominates")

        return "\n".join(patterns)

    def _calculate_nominal_mass(self, elements: dict) -> int:
        """Calculate nominal (integer) mass."""
        nominal_masses = {
            "H": 1, "C": 12, "N": 14, "O": 16, "F": 19, "Na": 23,
            "Si": 28, "P": 31, "S": 32, "Cl": 35, "Br": 79, "I": 127,
            "B": 11, "Al": 27, "Ca": 40, "Fe": 56,
        }
        total = 0
        for elem, count in elements.items():
            total += nominal_masses.get(elem, int(round(self._ISOTOPES.get(elem, [(elem, 1)])[0][0]))) * count
        return total

    def _check_nitrogen_rule(self, elements: dict, nominal_mass: int) -> str:
        """Apply the nitrogen rule."""
        n_count = elements.get("N", 0)
        is_odd = (nominal_mass % 2 == 1)

        if n_count == 0:
            if is_odd:
                return "⚠️ ODD mass with ZERO N atoms — check formula (radical cation exception possible)"
            else:
                return "✅ EVEN mass with zero N atoms — consistent with even-electron ion (or no N)"
        elif n_count % 2 == 1:  # odd number of N
            if is_odd:
                return "✅ ODD mass with ODD number of N atoms — consistent with nitrogen rule"
            else:
                return "⚠️ EVEN mass with ODD number of N — violates nitrogen rule (check formula)"
        else:  # even number of N
            if not is_odd:
                return "✅ EVEN mass with EVEN number of N atoms — consistent with nitrogen rule"
            else:
                return "⚠️ ODD mass with EVEN number of N — violates nitrogen rule (check formula)"

    def _run_base(self, formula: str) -> dict:
        """
        Calculate molecular ion data.

        Args:
            formula: Molecular formula or SMILES string

        Returns:
            Dict with complete molecular ion analysis
        """
        if not formula:
            raise ChemMCPError("Molecular formula or SMILES string is required.")

        parsed = self._parse_formula(formula)
        elements = parsed["elements"]

        if not elements:
            raise ChemMCPError(f"Could not parse formula: '{formula}'")

        # Calculate all properties
        exact_mass = self._calculate_exact_mass(elements)
        nominal_mass = self._calculate_nominal_mass(elements)
        isotope_pattern = self._calculate_isotope_pattern(elements)
        pattern_recognition = self._recognize_pattern(elements)
        nitrogen_rule = self._check_nitrogen_rule(elements, nominal_mass)

        # Build pretty formula string
        pretty_formula = "".join(
            f"{elem}{count if count > 1 else ''}"
            for elem, count in sorted(elements.items())
        )

        # Calculate Ring Double Bond Equivalent (RDBE / degree of unsaturation)
        rdbe = self._calculate_rdbe(elements)

        return {
            "ion_data": {
                "input": formula.strip(),
                "parsed_formula": pretty_formula,
                "element_counts": elements,
                "nominal_mass_mz": nominal_mass,
                "exact_mass": exact_mass,
                "mass_error_ppm": 0,  # Reference value
                "isotope_pattern": {
                    "M_percent": 100.0,
                    "M_plus_1_percent": isotope_pattern["M_plus_1_percent"],
                    "M_plus_2_percent": isotope_pattern["M_plus_2_percent"],
                    "M_plus_1_ratio": isotope_pattern["M_plus_1_relative_to_M"],
                    "M_plus_2_ratio": isotope_pattern["M_plus_2_relative_to_M"],
                },
                "isotope_pattern_recognition": pattern_recognition,
                "nitrogen_rule_check": nitrogen_rule,
                "rings_plus_double_bonds": rdbe,
                "unsaturation_interpretation": self._interpret_rdbe(rdbe),
                "notes": (
                    "Monoisotopic mass = sum of lightest isotope masses\n"
                    "M+1 mainly from ¹³C (~1.07% per C atom)\n"
                    "M+2 from ¹⁸O, ³⁴S, ³⁰Si, ³⁷Cl, ⁸¹Br\n"
                    "Nitrogen Rule: odd nominal mass → odd # of N atoms\n"
                    "RDBE = C - H/2 + N/2 + 1 (for CxHyNzOw)"
                ),
            }
        }

    def _calculate_rdbe(self, elements: dict) -> float:
        """Calculate Rings + Double Bonds Equivalent."""
        c = elements.get("C", 0)
        h = elements.get("H", 0)
        n = elements.get("N", 0)
        x = elements.get("X", 0)  # halogens treated as H
        # Count halogens as hydrogen equivalents
        halogens = sum(elements.get(x_elem, 0) for x_elem in ["Cl", "Br", "I", "F"])

        if c == 0:
            return 0.0

        rdbe = c - (h + halogens + n) / 2 + n / 2 + 1
        return round(rdbe, 1)

    def _interpret_rdbe(self, rdbe: float) -> str:
        """Interpret RDBE value."""
        if rdbe < 0:
            return f"⚠️ Negative RDBE ({rdbe}) — invalid formula or charged species!"
        elif rdbe == 0:
            return "Fully saturated acyclic compound (no rings, no double/triple bonds)"
        elif rdbe < 1:
            return f"One ring OR one π bond (RDBE = {rdbe})"
        else:
            int_part = int(rdbe)
            interpretations = []
            # Possible combinations
            rings = min(int_part, int_part)
            pi_bonds = int_part - rings
            interpretations.append(f"Total unsaturation: {int_part}")
            if int_part >= 4:
                interpretations.append(f"Possible: benzene ring (4 RDBE: 3 π bonds + 1 ring)")
            if int_part >= 2:
                interpretations.append(f"Possible: triple bond (2 RDBE) or double bond + ring or 2 double bonds")
            if int_part == 1:
                interpretations.append(f"Either one double bond OR one ring")
            return " | ".join(interpretations)

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        formula = input_params.strip()
        if not formula:
            raise ChemMCPError("Input required: molecular formula or SMILES")

        return self._run_base(formula)
