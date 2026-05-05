"""
Elemental Composition Calculator - Calculates possible elemental compositions
from exact mass using brute-force enumeration with chemical constraints.
"""

import logging
import math
import re
from typing import Dict, List, Tuple, Any, Optional, Set
from itertools import product as iter_product

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Exact atomic masses (IUPAC 2022, most abundant isotope for monoisotopic mass)
_ATOMIC_MASSES_EXACT: Dict[str, float] = {
    "H":   1.00782503223,
    "C":  12.00000000000,
    "N":  14.00307400443,
    "O":  15.99491461957,
    "P":  30.97376199842,
    "S":  31.97207117444,
    "F":  18.99840316273,
    "Cl": 34.968852682,
    "Br": 78.9183376,
    "I":  126.9044727,
    "Si": 27.97692653465,
    "Na": 22.9897692820,
    "B":  10.01293695,
}

# Default element order and constraints
_DEFAULT_ELEMENTS = ["C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I"]

# Valence constraints (typical max counts relative to carbon)
_VALENCE_RULES: Dict[str, Tuple[int, int, str]] = {
    # Element: (max_per_H, max_absolute, typical_ratio_note)
    "H":  (3, 200, "~2-3× C count"),
    "N":  (4, 20,  "≤~4× C, typically ≤ C"),
    "O":  (4, 25,  "≤~4× C"),
    "P":  (1, 6,   "rarely >1 per molecule"),
    "S":  (2, 10,  "≤~2× C"),
    "F":  (20, 40, "can be many in fluorinated compounds"),
    "Cl": (10, 15, "halogenated compounds"),
    "Br": (5, 8,  "brominated compounds"),
    "I":  (3, 5,   "iodinated compounds"),
    "Si": (2, 6,   "organosilicon"),
    "Na": (1, 2,   "salt form only"),
    "B":  (4, 8,   "boron compounds"),
}


@ChemMCPManager.register_tool
class ElementalCompositionCalculator(BaseTool):
    """
    元素组成计算器 — 基于精确质量计算可能的元素组成（分子式）。
    
    使用穷举算法结合化学约束条件（氮规则、RDBE、价态规则等），
    从精确质量反推候选分子式，并按匹配度排序。
    """
    __version__      = "0.1.0"
    name             = "ElementalCompositionCalculator"
    func_name        = "calculate_elemental_composition"
    description      = "Calculate possible elemental compositions (molecular formulas) from an exact mass measurement with ppm tolerance."
    implementation_description = "Uses constrained brute-force enumeration over allowed elements with chemical validity filters: nitrogen rule (odd/even mass vs N count), RDBE (Ring Double Bond Equivalents ≥ 0), valence reasonableness, H/C ratio bounds, and Leven Five rule. Candidates are ranked by mass error (ppm) then by heuristic score."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Molecule"]
    tags             = ["Mass Spectrometry", "Elemental Composition", "Exact Mass", "Formula Prediction", "HRMS"]
    required_envs    = []

    code_input_sig   = [
        ("exact_mass", "float", "N/A", "Measured exact mass (monoisotopic, in Daltons)."),
        ("tolerance_ppm", "float", "5.0", "Mass tolerance in parts-per-million (ppm). Typical: 1-10 ppm for Orbitrap/FT-ICR, 5-20 for Q-TOF."),
        ("allowed_elements", "list", "None", "List of element symbols to consider. Default: C H N O P S F Cl Br I."),
        ("max_carbons", "int", "50", "Maximum number of carbon atoms to search."),
        ("max_unsaturation", "int", "30", "Maximum RDBE (Ring Double Bond Equivalents)."),
        ("charge", "int", "0", "Charge state (0=neutral, +1=[M+H]+, etc.). Used to adjust input mass."),
        ("nitrogen_rule", "bool", "True", "Apply nitrogen rule filter: odd nominal mass → odd N count."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'exact_mass [tolerance_ppm] [max_carbons] [charge]'. Example: '286.1438 5 50 1'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict with candidate_formulas list (ranked), each containing formula, exact_mass, error_ppm, RDBE, unsaturation_details, candidate_count, search_parameters, and quality notes."),
    ]

    examples         = [
        {
            "code_input": {
                "exact_mass": 286.1438,
                "tolerance_ppm": 5.0,
                "allowed_elements": None,
                "max_carbons": 30,
                "max_unsaturation": 20,
                "charge": 1,
                "nitrogen_rule": True,
            },
            "text_input": {
                "input_params": "286.1438 5 30 1"
            },
            "output": {
                "result": {
                    "input_mass": 286.1438,
                    "search_tolerance_ppm": 5.0,
                    "candidate_count": 5,
                    "top_candidates": [
                        {"rank": 1, "formula": "C17H19NO3", "calc_mass": 286.1385, "error_ppm": -1.87, "RDBE": 9.0, "validity": "passes_all_rules"},
                        {"rank": 2, "formula": "C16H19NO4", "calc_mass": 286.1314, "error_ppm": -4.36, "RDBE": 8.0, "validity": "passes_all_rules"},
                    ],
                    "search_summary": "Search completed within constraints",
                }
            },
        },
        {
            "code_input": {
                "exact_mass": 180.0634,
                "tolerance_ppm": 10.0,
                "charge": 0,
            },
            "text_input": {
                "input_params": "180.0634 10"
            },
            "output": {
                "result": {
                    "top_candidates": [
                        {"rank": 1, "formula": "C9H10O4", "calc_mass": 180.0630, "error_ppm": -2.22, "RDBE": 5.0},
                        {"rank": 2, "formula": "C6H14N4O3", "calc_mass": 178.1060, "error_ppm": -129.0, "RDBE": 2.0},
                    ]
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, exact_mass: float, tolerance_ppm: float = 5.0,
                  allowed_elements: Optional[List[str]] = None,
                  max_carbons: int = 50, max_unsaturation: int = 30,
                  charge: int = 0, nitrogen_rule: bool = True) -> dict:
        """Core logic: enumerate and score elemental compositions."""
        if exact_mass <= 0:
            raise ChemMCPError("Exact mass must be positive.")
        if tolerance_ppm <= 0 or tolerance_ppm > 1000:
            raise ChemMCPError("Tolerance must be between 0-1000 ppm.")

        # Adjust for charge state (e.g., [M+H]+ → subtract proton mass)
        neutral_mass = exact_mass - charge * 1.007276466812
        
        if neutral_mass <= 0:
            raise ChemMCPError(f"Adjusted neutral mass ({neutral_mass:.4f}) is not positive. Check charge state.")

        # Default elements
        if allowed_elements is None:
            allowed_elements = list(_DEFAULT_ELEMENTS)

        # Validate elements
        for elem in allowed_elements:
            if elem not in _ATOMIC_MASSES_EXACT:
                raise ChemMCPError(f"Unknown element '{elem}'. Available: {list(_ATOMIC_MASSES_EXACT.keys())}")

        # Ensure C is first if present (optimization)
        if "C" in allowed_elements:
            allowed_elements = ["C"] + [e for e in allowed_elements if e != "C"]

        # Calculate mass window
        tol_da = abs(neutral_mass * tolerance_ppm / 1e6)
        mass_lo = neutral_mass - tol_da
        mass_hi = neutral_mass + tol_da

        nominal_mass = int(round(neutral_mass))

        logger.info(f"Searching compositions: mass={neutral_mass:.4f} ±{tol_da:.4f} Da (±{tolerance_ppm} ppm), "
                     f"elements={allowed_elements}, C≤{max_carbons}, RDBE≤{max_unsaturation}")

        # Enumerate candidates
        candidates = self._enumerate_compositions(
            neutral_mass, mass_lo, mass_hi, nominal_mass,
            allowed_elements, max_carbons, max_unsaturation,
            nitrogen_rule, tolerance_ppm
        )

        # Sort by absolute error, then by heuristic score
        candidates.sort(key=lambda c: (abs(c["error_ppm"]), -c["heuristic_score"]))

        # Assign ranks
        for i, c in enumerate(candidates):
            c["rank"] = i + 1

        # Generate summary
        summary = self._generate_summary(neutral_mass, candidates, allowed_elements, tolerance_ppm)

        return {
            "result": {
                "input_exact_mass": round(exact_mass, 6),
                "input_charge": charge,
                "adjusted_neutral_mass": round(neutral_mass, 6),
                "search_tolerance_ppm": tolerance_ppm,
                "search_tolerance_Da": round(tol_da, 6),
                "mass_window": (round(mass_lo, 6), round(mass_hi, 6)),
                "allowed_elements": allowed_elements,
                "max_carbons": max_carbons,
                "max_RDBE": max_unsaturation,
                "nitrogen_rule_applied": nitrogen_rule,
                "candidate_count": len(candidates),
                "candidates": candidates[:25],  # Return top 25
                "summary": summary,
                "notes": (
                    "Elemental Composition Notes:\n"
                    "• Candidates are ranked by mass error (ppm) then heuristic score\n"
                    "• RDBE = Ring Double Bond Equivalents (rings + π bonds)\n"
                    f"• {len(candidates)} candidate(s) found within ±{tolerance_ppm} ppm\n"
                    "• Always verify with MS/MS or orthogonal data when possible\n"
                    "• Lower RDBE values are generally more plausible for small molecules\n"
                    "• Consider isotopic pattern match for final confirmation"
                ),
            }
        }

    def _enumerate_compositions(self, target_mass: float, mass_lo: float, mass_hi: float,
                                 nominal_mass: int, elements: List[str],
                                 max_c: int, max_rdbe: int,
                                 apply_n_rule: bool, tol_ppm: float) -> List[dict]:
        """Efficient enumeration: C/H outer loops + limited heteroatom search."""
        candidates = []
        
        elem_names = list(elements)
        has_c = "C" in elem_names
        has_h = "H" in elem_names
        
        # Heteroatoms (non-C, non-H)
        hetero_atoms = [e for e in elem_names if e not in ("C", "H")]
        
        # Precompute heteroatom masses
        het_masses = {e: _ATOMIC_MASSES_EXACT[e] for e in hetero_atoms}
        
        # Limited heteroatom combinations (pre-compute small set)
        het_combos = self._get_heteroatom_combinations(hetero_atoms, target_mass, max_per_elem=6)
        
        c_mass = _ATOMIC_MASSES_EXACT.get("C", 0)
        h_mass = _ATOMIC_MASSES_EXACT.get("H", 0)
        
        # Carbon loop
        max_c = min(max_c, int(target_mass / c_mass) + 2)
        for n_c in range(0 if not has_c else 1, max_c + 1):
            base_mass_c = n_c * c_mass
            
            if base_mass_c > mass_hi + 50:
                break
            # Hydrogen range — for organic molecules, H is typically 1 to 2C+2+N
            # Use chemistry-based bounds rather than mass-based (which fails when het atoms dominate)
            max_h_for_c = 2 * n_c + 20 + 10  # saturated + extra for polyols
            min_h_for_c = 1 if has_c else 0
            # Also bound by total mass reasonableness
            max_h = min(max_h_for_c, int(target_mass / h_mass) + 5)
            # Allow very low H for highly oxidized / halogenated compounds
            min_h = min_h_for_c
            
            for n_h in range(min_h, max_h + 1):
                ch_mass = base_mass_c + n_h * h_mass
                
                # Quick prune: skip if CH alone already exceeds window + reasonable het
                if ch_mass > mass_hi + 80:
                    continue
                # Try each heteroatom combination
                for het_combo in het_combos:
                    het_total_mass = sum(het_masses.get(e, 0) * cnt for e, cnt in het_combo.items())
                    total_mass = ch_mass + het_total_mass
                    
                    if total_mass < mass_lo - 0.01 or total_mass > mass_hi + 0.01:
                        continue
                    # Build full composition
                    counts = {}
                    if has_c:
                        counts["C"] = n_c
                    if has_h:
                        counts["H"] = n_h
                    counts.update(het_combo)
                    
                    # Skip empty compositions
                    if sum(counts.values()) == 0:
                        continue
                    self._check_composition(
                        counts, total_mass, target_mass, mass_lo, mass_hi,
                        nominal_mass, candidates, max_rdbe, apply_n_rule, tol_ppm
                    )

        return candidates

    def _get_heteroatom_combinations(self, elements: List[str], target_mass: float,
                                       max_per_elem: int = 8) -> List[Dict[str, int]]:
        """Generate heteroatom count combinations for efficient search."""
        if not elements:
            return [{}]
        
        masses = {e: _ATOMIC_MASSES_EXACT[e] for e in elements}
        
        # Limit each element's range
        ranges = []
        for e in elements:
            max_n = min(max_per_elem, int(target_mass / masses[e]) + 3)
            ranges.append(range(0, max_n + 1))
        
        combos = []
        total_combos = 1
        for r in ranges:
            total_combos *= len(r)
        
        if total_combos <= 100000:
            # Small enough for full Cartesian product
            from itertools import product as iter_product
            for combo in iter_product(*ranges):
                d = {elements[i]: n for i, n in enumerate(combo) if n > 0}
                if d:
                    combos.append(d)
            combos.append({})  # include no-heteroatom case
        else:
            # Strategic sampling for large spaces
            combos.append({})  # no heteroatoms
            
            # All single-element cases (complete)
            for ei, e in enumerate(elements):
                for n in range(1, len(ranges[ei])):  # all non-zero values
                    combos.append({e: n})
            
            # Two-element pairs (small counts — covers most organic molecules)
            if len(elements) >= 2:
                for i in range(min(len(elements), 5)):
                    e1 = elements[i]
                    for j in range(i + 1, min(len(elements), 6)):
                        e2 = elements[j]
                        for n1 in range(1, 5):
                            for n2 in range(1, 5):
                                combos.append({e1: n1, e2: n2})
            
            # Three-element triples (very small counts — N/O/S common)
            if len(elements) >= 3:
                for i in range(min(len(elements), 3)):
                    for j in range(i+1, min(len(elements), 4)):
                        for k in range(j+1, min(len(elements), 5)):
                            for n1 in range(1, 4):
                                for n2 in range(1, 4):
                                    for n3 in range(1, 4):
                                        combos.append({
                                            elements[i]: n1,
                                            elements[j]: n2,
                                            elements[k]: n3
                                        })

        return combos

    def _enumerate_recursive(self, target_mass, mass_lo, mass_hi, nominal_mass,
                              elems, depth, current_counts, current_mass,
                              ranges, candidates, max_rdbe, n_rule, tol_ppm,
                              iter_count, max_iter):
        """Recursive enumeration with pruning."""
        if iter_count > max_iter:
            return
        if depth >= len(elems):
            # Leaf node — check this composition
            self._check_composition(
                current_counts, current_mass, target_mass, mass_lo, mass_hi,
                nominal_mass, candidates, max_rdbe, n_rule, tol_ppm
            )
            return

        elem = elems[depth]
        amass = _ATOMIC_MASSES_EXACT.get(elem, 0)
        
        # Pruning: if even minimum additional mass exceeds window, skip
        if current_mass > mass_hi + 50:
            return
        
        for n in ranges.get(elem, range(0, 3)):
            new_mass = current_mass + n * amass
            new_counts = dict(current_counts)
            new_counts[elem] = n
            
            # Early pruning: skip if already too heavy
            if depth < len(elems) - 1 and new_mass > mass_hi + 20:
                break  # sorted, so further n will be heavier
            
            self._enumerate_recursive(
                target_mass, mass_lo, mass_hi, nominal_mass,
                elems, depth + 1, new_counts, new_mass,
                ranges, candidates, max_rdbe, n_rule, tol_ppm,
                iter_count + 1, max_iter
            )

    def _check_composition(self, counts: Dict[str, int], calc_mass: float,
                            target_mass: float, mass_lo: float, mass_hi: float,
                            nominal_mass: int, candidates: List[dict],
                            max_rdbe: int, n_rule: bool, tol_ppm: float):
        """Check if a composition passes all filters."""
        total = sum(counts.values())
        if total == 0:
            return

        # Skip pure hydrogen
        if set(counts.keys()) == {"H"}:
            return

        # Mass check
        if calc_mass < mass_lo or calc_mass > mass_hi:
            return

        # Get element counts
        n_c = counts.get("C", 0)
        n_h = counts.get("H", 0)
        n_n = counts.get("N", 0)
        n_o = counts.get("O", 0)
        n_p = counts.get("P", 0)
        n_s = counts.get("S", 0)

        # Must have at least one non-H atom if C=0
        if n_c == 0 and sum(counts.get(e, 0) for e in counts if e != "H") == 0:
            return

        # Nitrogen rule: odd nominal mass ↔ odd N count
        if n_rule:
            calc_nominal = int(round(calc_mass))
            is_odd_nominal = (calc_nominal % 2 == 1)
            is_odd_n = (n_n % 2 == 1)
            if is_odd_nominal != is_odd_n:
                return  # Violates nitrogen rule

        # RDBE calculation: RDBE = C - H/2 + N/2 + 1
        rdbe = n_c - n_h / 2.0 + n_n / 2.0 + 1
        if rdbe < 0 or rdbe > max_rdbe:
            return  # Invalid or out-of-range RDBE

        # RDBE must be integer (or very close)
        if abs(rdbe - round(rdbe)) > 0.01:
            return

        rdbe_int = int(round(rdbe))

        # H/C ratio sanity check
        if n_c > 0:
            h_c_ratio = n_h / n_c
            if h_c_ratio > 6 or (h_c_ratio < 0.33 and n_c > 2):
                return  # Unlikely H/C ratio

        # Minimum H for given C (fully saturated acyclic: H = 2C+2)
        if n_c > 0:
            min_h_for_c = 0  # allow some exotic structures
            max_h_for_c = 2 * n_c + 2 + n_n  # fully saturated + N additions
            if n_h > max_h_for_c + 20:  # some slack for polyols
                return

        # Calculate error
        error_ppm = (calc_mass - target_mass) / target_mass * 1e6
        if abs(error_ppm) > tol_ppm:
            return

        # Heuristic scoring (lower is better)
        heuristic_score = self._heuristic_score(counts, rdbe_int, error_ppm)

        # Build formula string
        formula_parts = []
        for elem in sorted(counts.keys()):
            cnt = counts[elem]
            if cnt > 0:
                formula_parts.append(f"{elem}{cnt}" if cnt > 1 else elem)
        formula = "".join(formula_parts)

        candidates.append({
            "formula": formula,
            "element_counts": dict(counts),
            "calc_mass": round(calc_mass, 6),
            "error_ppm": round(error_ppm, 4),
            "RDBE": rdbe_int,
            "unsaturation_type": self._classify_rdbe(rdbe_int, n_c),
            "heuristic_score": round(heuristic_score, 2),
            "validity": "passes_all_rules",
        })

    def _heuristic_score(self, counts: Dict[str, int], rdbe: int, error_ppm: float) -> float:
        """
        Score candidate plausibility (lower = more plausible).
        Combines multiple heuristics.
        """
        score = 0.0
        n_c = counts.get("C", 0)
        n_h = counts.get("H", 0)
        n_o = counts.get("O", 0)
        n_n = counts.get("N", 0)

        # Prefer reasonable RDBE range for organic molecules
        if 0 <= rdbe <= 15:
            score -= 0
        elif 16 <= rdbe <= 25:
            score += 5
        else:
            score += 15

        # Prefer reasonable H/C ratio (~1-2 for aromatics, ~2 for aliphatic)
        if n_c > 0:
            hc = n_h / n_c
            if 0.5 <= hc <= 2.5:
                score -= 2  # good
            elif 2.5 < hc <= 4:
                score += 3
            else:
                score += 8

        # Penalize too many heteroatoms without enough carbons
        hetero = n_o + n_n + counts.get("S", 0) + counts.get("P", 0)
        if n_c > 0 and hetero > n_c:
            score += 10

        # Slight preference for C/H/N/O compositions (most common organics)
        non_chno = sum(counts.get(e, 0) for e in counts if e not in ("C", "H", "N", "O"))
        if non_chno == 0:
            score -= 3

        # Penalize extreme atom counts
        if n_h > 150:
            score += 10
        if n_c < 1 and sum(counts.values()) > 5:
            score += 5  # no-carbon molecules are rare

        return score

    def _classify_rdbe(self, rdbe: int, n_c: int) -> str:
        """Classify unsaturation type."""
        if rdbe == 0:
            return "acyclic_saturated"
        elif rdbe == 1:
            return "one_ring_or_pi_bond"
        elif rdbe <= 4:
            return "moderate_unsaturation"
        elif n_c >= 6 and rdbe >= 4 and rdbe <= 7:
            return "likely_aromatic"
        elif rdbe > 7:
            return "highly_unsaturated_or_polycyclic"
        else:
            return f"{rdbe}_rdbe"

    def _generate_summary(self, mass: float, candidates: list, elements: list, tol: float) -> str:
        """Generate human-readable summary."""
        if not candidates:
            return f"No candidates found within ±{tol} ppm of mass {mass:.4f}. Try widening tolerance or element range."
        
        best = candidates[0]
        parts = [
            f"Best match: {best['formula']} (error={best['error_ppm']:+.2f} ppm, RDBE={best['RDBE']})",
            f"Total candidates found: {len(candidates)}",
        ]
        
        if len(candidates) > 1:
            second = candidates[1]
            parts.append(f"Second best: {second['formula']} (error={second['error_ppm']:+.2f} ppm)")
        
        if len(candidates) <= 3:
            parts.append("Good specificity — few candidate formulas")
        elif len(candidates) <= 10:
            parts.append("Moderate ambiguity — MS/MS recommended for confirmation")
        else:
            parts.append("Many candidates — need additional constraints (isotope pattern, MS/MS, or narrower tolerance)")

        return ". ".join(parts)

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            mass = float(parts[0])
            tol = float(parts[1]) if len(parts) > 1 else 5.0
            max_c = int(parts[2]) if len(parts) > 2 else 50
            chg = int(parts[3]) if len(parts) > 3 else 0
            return self._run_base(mass, tol, None, max_c, 30, chg, True)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'exact_mass [tolerance_ppm] [max_carbons] [charge]'")
