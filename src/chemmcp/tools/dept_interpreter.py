"""
DEPT Interpreter - interprets DEPT (Distortionless Enhancement by Polarization Transfer)
NMR spectra to distinguish CH3, CH2, CH, and quaternary carbons.
"""

import logging
from typing import Dict, List, Tuple, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class DeptInterpreter(BaseTool):
    __version__      = "0.1.0"
    name             = "DeptInterpreter"
    func_name        = "interpret_dept"
    description      = "Interpret DEPT NMR spectra to distinguish CH₃ (positive), CH₂ (negative), CH (positive), and quaternary Cq (absent in DEPT) carbons."
    implementation_description = "Compares DEPT-90 and DEPT-135 spectra against regular 13C NMR to classify each carbon signal as CH3/CH2/CH/Cq based on phase behavior."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Molecule"]
    tags             = ["NMR", "DEPT", "Spectroscopy", "Carbon-13"]
    required_envs    = []

    code_input_sig   = [
        ("dept_90_peaks", "list", "[]", "List of chemical shifts (ppm) observed in DEPT-90 spectrum (only CH signals appear as positive peaks)."),
        ("dept_135_peaks", "list", "[]", "List of chemical shifts (ppm) observed in DEPT-135 spectrum (CH & CH₃ positive, CH₂ negative)."),
        ("regular_13c_peaks", "list", "N/A", "Complete list of all ¹³C NMR chemical shifts (ppm). Required for identifying quaternary carbons."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Semicolon-separated: 'dept90_shifts;dept135_shifts;all_13c_shifts'. Example: '20,30,70;10(+),15(+),20(-),30(-),40(+),50(+),70(+);10,15,20,25,30,40,50,60,70'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict containing: carbon_assignments (list per carbon with shift, type, rationale), summary_counts (CH3/CH2/CH/Cq counts), interpretation_notes."),
    ]

    examples         = [
        {
            "code_input": {
                "dept_90_peaks": [28.0, 68.0],
                "dept_135_peaks": [14.0, 22.5, 28.0, 31.5, 42.0, 68.0],
                "regular_13c_peaks": [14.0, 22.5, 28.0, 31.5, 42.0, 60.0, 68.0, 210.0],
            },
            "text_input": {"input_params": "28,68;14,22.5,28,31.5,42,68;14,22.5,28,31.5,42,60,68,210"},
            "output": {
                "result": {
                    "carbon_assignments": [
                        {"shift": 14.0, "type": "CH3", "rationale": "In DEPT-135 only (not DEPT-90): positive → CH₃"},
                        {"shift": 22.5, "type": "CH2", "rationale": "In DEPT-135 only (not DEPT-90): negative → CH₂"},
                        {"shift": 28.0, "type": "CH", "rationale": "In both DEPT-90 and DEPT-135: positive → CH"},
                        {"shift": 31.5, "type": "CH2", "rationale": "In DEPT-135 only (not DEPT-90): negative → CH₂"},
                        {"shift": 42.0, "type": "CH", "rationale": "In both DEPT-90 and DEPT-135: positive → CH"},
                        {"shift": 60.0, "type": "Cq", "rationale": "Absent from both DEPT-90 and DEPT-135 but present in ¹³C → Quaternary"},
                        {"shift": 68.0, "type": "CH", "rationale": "In both DEPT-90 and DEPT-135: positive → CH"},
                        {"shift": 210.0, "type": "Cq", "rationale": "Absent from both DEPT-90 and DEPT-135 but present in ¹³C → Quaternary (likely C=O)"},
                    ],
                    "summary_counts": {"CH3": 1, "CH2": 2, "CH": 3, "Cq": 2},
                    "interpretation_notes": ["C=O at 210 ppm confirmed as quaternary", "O-bearing CH at 68 ppm"],
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, dept_90_peaks: list, dept_135_peaks: list, regular_13c_peaks: list) -> dict:
        """Core logic: classify each 13C peak using DEPT information."""
        if not regular_13c_peaks:
            raise ChemMCPError("Regular 13C NMR peak list is required.")

        # Normalize inputs: handle both plain lists and signed lists for DEPT-135
        dept_90_set = set(round(float(p), 3) for p in dept_90_peaks)
        dept_135_positive = set()
        dept_135_negative = set()

        for p in dept_135_peaks:
            ps = str(p).strip()
            # Check for sign notation like "(+)" or "(-)"
            if "(" in ps or ps.startswith("+") or ps.startswith("-"):
                import re
                m = re.match(r'([+-]?\d+\.?\d*)\s*[\(\[]([+-])[\)\]]', ps)
                if m:
                    val = round(float(m.group(1)), 3)
                    sign = m.group(2)
                    if sign == '+':
                        dept_135_positive.add(val)
                    else:
                        dept_135_negative.add(val)
                else:
                    val = round(float(ps.replace('(', '').replace(')', '').replace('+', '').replace('-', '').strip()), 3)
                    dept_135_positive.add(val)
            else:
                val = round(float(p), 3)
                dept_135_positive.add(val)

        all_13c = [round(float(p), 3) for p in regular_13c_peaks]

        assignments = []
        notes = []

        for shift in sorted(all_13c):
            in_90 = shift in dept_90_set
            in_135_pos = shift in dept_135_positive
            in_135_neg = shift in dept_135_negative

            if in_90:
                ctype = "CH"
                rationale = "In both DEPT-90 and DEPT-135: positive → CH"
            elif in_135_pos:
                ctype = "CH3"
                rationale = "In DEPT-135 only (not DEPT-90): positive → CH₃"
            elif in_135_neg:
                ctype = "CH2"
                rationale = "In DEPT-135 only (not DEPT-90): negative → CH₂"
            elif not in_90 and shift not in dept_135_positive and shift not in dept_135_negative:
                ctype = "Cq"
                rationale = "Absent from both DEPT-90 and DEPT-135 but present in ¹³C → Quaternary"
            else:
                ctype = "Unknown"
                rationale = "Unexpected DEPT pattern"

            assignments.append({
                "shift": shift,
                "type": ctype,
                "rationale": rationale,
            })

            # Add interpretive notes for characteristic regions
            if ctype == "Cq" and shift > 160:
                notes.append(f"Quaternary at {shift} ppm likely a carbonyl (C=O)")
            elif ctype == "Cq" and 110 <= shift <= 160:
                notes.append(f"Quaternary at {shift} ppm may be an aromatic or alkene C without attached H")
            elif ctype == "CH" and 55 <= shift <= 90:
                notes.append(f"CH at {shift} ppm suggests C bonded to electronegative atom (O/N/halogen)")
            elif ctype == "CH3" and shift < 15:
                notes.append(f"CH₃ at {shift} ppm typical of terminal methyl group")

        # Summary counts
        counts: Dict[str, int] = {"CH3": 0, "CH2": 0, "CH": 0, "Cq": 0}
        for a in assignments:
            t = a["type"]
            if t in counts:
                counts[t] += 1

        return {
            "carbon_assignments": assignments,
            "summary_counts": counts,
            "interpretation_notes": notes,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split(";")

            def parse_peak_list(s):
                s = s.strip()
                if not s:
                    return []
                items = [x.strip() for x in s.split(",")]
                result = []
                for item in items:
                    if item:
                        result.append(item)
                return result

            dept_90 = parse_peak_list(parts[0]) if len(parts) > 0 else []
            dept_135 = parse_peak_list(parts[1]) if len(parts) > 1 else []
            all_13c = parse_peak_list(parts[2]) if len(parts) > 2 else []

            return self._run_base(dept_90, dept_135, all_13c)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'dept90_shifts;dept135_shifts;all_13c_shifts'")
