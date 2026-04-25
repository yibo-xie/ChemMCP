"""
Spectrum to Structure Inference - guided analysis combining IR, NMR, and MS data
to infer molecular structure features.
"""

import logging
from typing import Dict, List, Tuple, Any, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# IR peak reference: wavenumber (cm-1) -> functional group interpretation
IR_REFERENCE: List[Dict[str, Any]] = [
    {"range": (3650, 3200), "strength": "broad", "fg": "O-H stretch (alcohol/phenol)", "notes": "Broad, H-bonded"},
    {"range": (3500, 3300), "strength": "medium", "fg": "N-H stretch (amine/amide)", "notes": "1°=2 peaks, 2°=1 peak"},
    {"range": (3300, 3000), "strength": "medium", "fg": "≡C-H stretch (alkyne)", "notes": "Sharp (~3300 cm⁻¹)"},
    {"range": (3100, 3000), "strength": "medium", "fg": "=C-H stretch (alkene/aromatic)", "notes": ">3000 cm⁻¹"},
    {"range": (3000, 2850), "strength": "strong", "fg": "C-H stretch (alkane)", "notes": "<3000 cm⁻¹"},
    {"range": (2750, 2650), "strength": "weak", "fg": "C-H stretch (aldehyde)", "notes": "Doublet (Fermi resonance)"},
    {"range": (2250, 2220), "strength": "medium", "fg": "C≡N stretch (nitrile)", "notes": "Sharp peak"},
    {"range": (2260, 2190), "strength": "weak-medium", "fg": "C≡C stretch (alkyne)", "notes": "Terminal alkyne stronger"},
    {"range": (1760, 1690), "strength": "strong", "fg": "C=O stretch (carbonyl)", "notes": "Very characteristic"},
    {"range": (1680, 1640), "strength": "variable", "fg": "C=C stretch (alkene)", "notes": "Medium intensity"},
    {"range": (1620, 1600), "strength": "variable", "fg": "C=C aromatic ring", "notes": "Often with ~1500 cm⁻¹ companion"},
    {"range": (1510, 1450), "strength": "variable", "fg": "C=C aromatic ring / C-H bend", "notes": "Fingerprint region start"},
    {"range": (1465, 1440), "strength": "medium", "fg": "C-H scissoring (CH₂)", "notes": ""},
    {"range": (1380, 1370), "strength": "medium", "fg": "C-H umbrella (CH₃)", "notes": ""},
    {"range": (1300, 1000), "strength": "strong", "fg": "C-O stretch (alcohol/ether/ester)", "notes": "Strong, broad region"},
    {"range": (980, 850), "strength": "strong", "fg": "=C-H bend (alkene)", "notes": "Trans ~965 cm⁻¹"},
    {"range": (770, 710), "strength": "strong", "fg": "C-H bend (aromatic out-of-plane)", "notes": "Substitution pattern indicator"},
]

# Carbonyl sub-region details
CARBONYL_DETAIL: Dict[str, Dict[str, Any]] = {
    "saturated_ester":     {"range": (1750, 1735), "typical": 1740},
    "unsaturated_ester":   {"range": (1730, 1715), "typical": 1722},
    "aldehyde":            {"range": (1740, 1720), "typical": 1730},
    "ketone":              {"range": (1725, 1705), "typical": 1715},
    "carboxylic_acid":     {"range": (1725, 1700), "typical": 1710},
    "conjugated_ketone":   {"range": (1700, 1680), "typical": 1685},
    "amide":               {"range": (1690, 1630), "typical": 1660},
    "alpha_beta_unsat":    {"range": (1705, 1680), "typical": 1693},
}

# NMR chemical shift ranges (ppm) -> proton type
H_NMR_RANGES: List[Dict[str, Any]] = [
    {"range": (9.5, 10.0), "type": "Aldehyde H", "multiplicity_common": "s"},
    {"range": (6.5, 8.5),  "type": "Aromatic H", "multiplicity_common": "m"},
    {"range": (5.5, 6.5),  "type": "Vinylic H (=C-H)", "multiplicity_common": "m"},
    {"range": (4.5, 5.5),  "type": "Alkynyl H or O/N attached C-H", "multiplicity_common": "s/t"},
    {"range": (3.3, 4.5),  "type": "H on C bonded to electronegative atom (O, N, halogen)", "multiplicity_common": "q/s/t/m"},
    {"range": (2.3, 3.0),  "type": "α to carbonyl / benzylic / propargylic", "multiplicity_common": "q/q'"},
    {"range": (2.1, 2.5),  "type": "Acetyl CH₃ (COCH₃)", "multiplicity_common": "s"},
    {"range": (1.8, 2.3),  "type": "Allylic / α to unsaturation / benzylic CH₂", "multiplicity_common": "m/t"},
    {"range": (1.2, 1.6),  "type": "Aliphatic CH₂", "multiplicity_common": "m/pent"},
    {"range": (0.85, 1.15), "type": "Aliphatic CH₃", "multiplicity_common": "t/d"},
]

# C-13 NMR shift ranges
C13_NMR_RANGES: List[Dict[str, Any]] = [
    {"range": (165, 220), "type": "Carbonyl C (ketone/aldehyde)"},
    {"range": (155, 180), "type": "Carboxylic acid / ester / amide carbonyl"},
    {"range": (120, 150), "type": "Aromatic / alkene C"},
    {"range": (100, 130), "type": "Alkene C"},
    {"range": (60, 90),   "type": "C bonded to O (alcohol, ether)"},
    {"range": (40, 65),   "type": "C-N / C-halogen"},
    {"range": (25, 50),   "type": "Aliphatic C-C / α to functional group"},
    {"range": (10, 25),   "type": "Aliphatic CH₃ / CH₂ (remote from FG)"},
]


@ChemMCPManager.register_tool
class SpectrumToStructure(BaseTool):
    __version__      = "0.1.0"
    name             = "SpectrumToStructure"
    func_name        = "analyze_spectrum_to_structure"
    description      = "Guided analysis to infer molecular structure from combined spectral data (IR + ¹H NMR + MS)."
    implementation_description = "Uses rule-based pattern matching against IR wavenumber databases, NMR chemical shift ranges, and molecular formula analysis to suggest structural features."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Molecule"]
    tags             = ["Spectroscopy", "Structure Elucidation", "IR", "NMR", "MS"]
    required_envs    = []

    code_input_sig   = [
        ("molecular_formula", "str", "N/A", "Molecular formula (e.g., 'C6H12O')."),
        ("ir_peaks", "list", "[]", "List of IR peaks as [(wavenumber, intensity), ...]. Intensity: 'strong', 'medium', 'weak'."),
        ("nmr_peaks", "list", "[]", "List of ¹H NMR peaks as [(chemical_shift_ppm, multiplicity, integration), ...]. Multiplicity: s/d/t/q/dd/m/etc."),
        ("ms_mz", "float", "0", "Molecular ion M+ peak m/z value (0 if not provided)."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Semicolon-separated: 'formula;IR(wavenumber:intensity,...);NMR(shift:mult:integration,...);MS_mz'. Example: 'C6H12O;(1715,strong);(2.1,s,3H);102'"),
    ]

    output_sig       = [
        ("result", "dict", "Comprehensive analysis report including: molecular_info, ir_analysis, nmr_analysis, ms_analysis, degree_of_unsaturation, suggested_functional_groups, structural_hints, next_steps."),
    ]

    examples         = [
        {
            "code_input": {
                "molecular_formula": "C6H12O",
                "ir_peaks": [[1715, "strong"], [2950, "strong"]],
                "nmr_peaks": [[2.1, "s", 3], [2.4, "t", 2], [1.6, "sextet", 2], [0.95, "t", 3]],
                "ms_mz": 100.0,
            },
            "text_input": {"input_params": "C6H12O;(1715,strong),(2950,strong);(2.1,s,3),(2.4,t,2),(1.6,sextet,2),(0.95,t,3);100"},
            "output": {
                "result": {
                    "molecular_formula": "C6H12O",
                    "degree_of_unsaturation": 1,
                    "suggested_functional_groups": ["Ketone (C=O)", "Alkane chains"],
                    "structural_hints": ["Likely a ketone: one degree of unsaturation + strong IR at 1715 cm⁻¹", "Ethyl groups suggested by triplet-quartet patterns"],
                    "next_steps": ["Confirm with 2D NMR (COSY, HSQC) for connectivity", "Compare with known spectra in database"],
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    @staticmethod
    def _calc_dou(formula: str) -> int:
        """Calculate degrees of unsaturation (hydrogen deficiency index)."""
        import re
        elem_counts: Dict[str, int] = {}
        for elem, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
            if elem:
                elem_counts[elem] = elem_counts.get(elem, 0) + (int(count) if count else 1)

        C = elem_counts.get("C", 0)
        H = elem_counts.get("H", 0)
        N = elem_counts.get("N", 0)
        X = elem_counts.get("Cl", 0) + elem_counts.get("Br", 0) + elem_counts.get("F", 0)
        # Each halogen counts as 1 H for DOU calculation
        dou = C + 1 - (H + X - N) / 2
        return max(0, int(dou))

    def _run_base(self, molecular_formula: str, ir_peaks: list = None, nmr_peaks: list = None, ms_mz: float = 0) -> dict:
        """Core logic: combined spectral analysis."""
        if not molecular_formula:
            raise ChemMCPError("Molecular formula is required.")

        ir_peaks = ir_peaks or []
        nmr_peaks = nmr_peaks or []

        result: Dict[str, Any] = {
            "molecular_formula": molecular_formula,
            "degree_of_unsaturation": self._calc_dou(molecular_formula),
            "ir_analysis": self._analyze_ir(ir_peaks),
            "nmr_analysis": self._analyze_nmr(nmr_peaks),
            "ms_analysis": self._analyze_ms(ms_mz, molecular_formula),
            "suggested_functional_groups": [],
            "structural_hints": [],
            "next_steps": [],
        }

        # Synthesize findings
        fg_list = result["ir_analysis"].get("identified_groups", [])
        fg_list += result["nmr_analysis"].get("inferred_proton_types", [])
        result["suggested_functional_groups"] = list(set(fg_list))

        dou = result["degree_of_unsaturation"]
        hints = []
        has_carbonyl = any("carbonyl" in g.lower() or "C=O" in g for g in fg_list)
        has_aromatic = any("aromatic" in g.lower() for g in fg_list)
        has_alkene = any("vinylic" in g.lower() or "alkene" in g.lower() for g in fg_list)

        if dou == 0 and not has_carbonyl:
            hints.append(f"DOU={dou}: saturated molecule (acyclic or cyclic without π bonds)")
        elif dou >= 4 and has_aromatic:
            hints.append(f"DOU={dou}: consistent with an aromatic ring (DOU ≥ 4)")
        elif has_carbonyl:
            remaining_dou = dou - 1
            if remaining_dou > 0:
                hints.append(f"DOU={dou}: one C=O accounts for 1 DOU, {remaining_dou} remaining (possibly rings or double bonds)")
            else:
                hints.append(f"DOU={dou}: fully accounted for by the carbonyl group")
        elif dou > 0:
            hints.append(f"DOU={dou}: indicates {dou} ring(s) and/or π bond(s)")

        # Check integration consistency
        total_H_from_nmr = sum(p[2] for p in nmr_peaks) if nmr_peaks else 0
        import re
        h_in_formula = 0
        for elem, cnt in re.findall(r"([A-Z][a-z]?)(\d*)", molecular_formula):
            if elem == "H":
                h_in_formula = int(cnt) if cnt else 1

        if total_H_from_nmr > 0 and h_in_formula > 0:
            ratio = total_H_from_nmr / h_in_formula
            if abs(ratio - 1.0) < 0.05:
                hints.append("✓ NMR integration matches formula hydrogen count")
            elif ratio < 1.0:
                hints.append(f"⚠ NMR integration ({total_H_from_nmr}H) < formula H ({h_in_formula}H): exchangeable protons (OH, NH) may not be visible")
            else:
                hints.append(f"⚠ NMR integration ({total_H_from_nmr}H) > formula H ({h_in_formula}H): check integration values")

        result["structural_hints"] = hints
        result["next_steps"] = [
            "Acquire DEPT-135 and 13C NMR for carbon skeleton confirmation",
            "Consider 2D experiments (COSY, HSQC, HMBC) for connectivity mapping",
            "Cross-reference with spectral database (SDBS, NIST)",
        ]

        return result

    def _analyze_ir(self, ir_peaks: list) -> dict:
        """Analyze IR peaks against reference database."""
        identified = []
        details = []
        has_carbonyl = False
        carbonyl_region = None

        for peak in ir_peaks:
            wn = float(peak[0])
            intensity = str(peak[1]) if len(peak) > 1 else "unknown"

            matched = False
            for ref in IR_REFERENCE:
                lo, hi = ref["range"]
                if lo <= wn <= hi:
                    identified.append(ref["fg"])
                    details.append({
                        "wavenumber": wn,
                        "intensity": intensity,
                        "assignment": ref["fg"],
                        "notes": ref["notes"],
                    })
                    matched = True
                    if "carbonyl" in ref["fg"].lower() or "C=O" in ref["fg"]:
                        has_carbonyl = True
                        carbonyl_region = wn
                    break

            if not matched:
                # Check carbonyl detail regions
                if 1690 <= wn <= 1760:
                    has_carbonyl = True
                    carbonyl_region = wn
                    best_match = None
                    for ctype, cdata in CARBONYL_DETAIL.items():
                        if cdata["range"][0] <= wn <= cdata["range"][1]:
                            best_match = ctype.replace("_", " ").title()
                            break
                    cname = best_match or f"Carbonyl (~{wn} cm⁻¹)"
                    identified.append(cname)
                    details.append({
                        "wavenumber": wn,
                        "intensity": intensity,
                        "assignment": cname,
                        "notes": "Carbonyl region - check sub-type",
                    })
                else:
                    details.append({
                        "wavenumber": wn,
                        "intensity": intensity,
                        "assignment": "Unassigned",
                        "notes": "No match in standard reference table",
                    })

        return {
            "peak_count": len(ir_peaks),
            "identified_groups": identified,
            "peak_details": details,
            "has_carbonyl": has_carbonyl,
            "carbonyl_wavenumber": carbonyl_region,
        }

    def _analyze_nmr(self, nmr_peaks: list) -> dict:
        """Analyze 1H NMR peaks."""
        inferred_types = []
        peak_details = []
        total_integration = 0

        for peak in nmr_peaks:
            shift = float(peak[0])
            mult = str(peak[1]) if len(peak) > 1 else "unknown"
            integ = float(peak[2]) if len(peak) > 2 else 0
            total_integration += integ

            matched_type = "Unknown"
            for ref in H_NMR_RANGES:
                lo, hi = ref["range"]
                if lo <= shift <= hi:
                    matched_type = ref["type"]
                    break

            inferred_types.append(matched_type)
            peak_details.append({
                "chemical_shift_ppm": round(shift, 2),
                "multiplicity": mult,
                "integration": integ,
                "inferred_type": matched_type,
            })

        return {
            "peak_count": len(nmr_peaks),
            "total_integration": total_integration,
            "inferred_proton_types": inferred_types,
            "peak_details": peak_details,
        }

    def _analyze_ms(self, ms_mz: float, formula: str) -> dict:
        """Basic MS analysis."""
        if ms_mz <= 0:
            return {"note": "No MS data provided"}

        import re
        # Calculate approximate MW from formula
        atomic_weights = {"C": 12.01, "H": 1.008, "O": 16.00, "N": 14.01, "S": 32.07, "P": 30.97}
        calc_mw = 0
        for elem, cnt in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
            mw = atomic_weights.get(elem, 0)
            calc_mw += mw * (int(cnt) if cnt else 1)

        nominal_ms = int(round(ms_mz))
        nominal_calc = int(round(calc_mw))

        return {
            "observed_m_z": ms_mz,
            "nominal_mass": nominal_ms,
            "calculated_mw_from_formula": round(calc_mw, 2),
            "mass_match": abs(nominal_ms - nominal_calc) <= 1,
            "nitrogen_rule_even": (nominal_ms % 2 == 0),
            "nitrogen_rule_note": "Even nominal mass suggests even number of N atoms (or zero)" if nominal_ms % 2 == 0 else "Odd nominal mass suggests odd number of N atoms",
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split(";")
            formula = parts[0].strip()

            ir_peaks = []
            nmr_peaks = []
            ms_mz = 0.0

            for part in parts[1:]:
                part = part.strip()
                if part.upper().startswith("IR") or part.startswith("(") and "," in part:
                    # Parse IR peaks like (1715,strong) or IR:(1715,strong)
                    inner = part
                    if ":" in inner:
                        inner = inner.split(":", 1)[1]
                    for m in re.findall(r'\(([^)]+)\)', inner):
                        items = m.split(",")
                        if len(items) >= 2:
                            ir_peaks.append([float(items[0].strip()), items[1].strip()])
                elif part.upper().startswith("NMR") or ":" in part and any(c.isdigit() for c in part.split(":")[1][:3]):
                    inner = part
                    if ":" in inner:
                        inner = inner.split(":", 1)[1]
                    for m in re.findall(r'\(([^)]+)\)', inner):
                        items = m.split(",")
                        if len(items) >= 3:
                            nmr_peaks.append([float(items[0].strip()), items[1].strip(), float(items[2].strip())])
                else:
                    try:
                        val = float(part.strip())
                        if val > 20:  # likely MS m/z
                            ms_mz = val
                    except ValueError:
                        pass

            return self._run_base(formula, ir_peaks, nmr_peaks, ms_mz)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'formula;IR(wn,int);NMR(shift,mult,int);mz'")
