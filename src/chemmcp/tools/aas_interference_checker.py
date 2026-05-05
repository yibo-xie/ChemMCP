"""
AAS Interference Checker — AAS光谱干扰和化学干扰诊断工具 (#325)

功能：
  诊断原子吸收光谱(AAS)测量中可能遇到的光谱干扰、化学干扰和电离干扰，
  并提供消除/缓解建议。
"""

import logging
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ── 光谱线重叠数据库 (analyte: [interferent, interferent_wavelength, separation_pm]) ──
# 数据来源: NIST Atomic Spectra Database, Dean's Analytical Chemistry Handbook
SPECTRAL_OVERLAP_DB: Dict[str, List[Dict[str, Any]]] = {
    "cd": [
        {"interferent": "Fe", "interferent_line_nm": 228.802, "analyte_line_nm": 228.802,
         "separation_pm": 0, "severity": "critical",
         "note": "Direct line overlap — Fe 228.802 overlaps Cd 228.802 exactly."},
    ],
    "zn": [
        {"interferent": "Fe", "interferent_line_nm": 213.859, "analyte_line_nm": 213.856,
         "separation_pm": 3, "severity": "high",
         "note": "Fe 213.859 is only 3 pm from Zn 213.856 — high-resolution or alternate line needed."},
    ],
    "al": [
        {"interferent": "CaOH", "interferent_line_nm": 309.3, "analyte_line_nm": 309.271,
         "separation_pm": 29, "severity": "medium",
         "note": "Molecular band interference from CaOH in nitrous oxide flame."},
    ],
    "hg": [
        {"interferent": "Fe", "interferent_line_nm": 253.682, "analyte_line_nm": 253.652,
         "separation_pm": 30, "severity": "medium",
         "note": "Fe nearby line; use cold vapor technique to avoid."},
        {"interferent": "Co", "interferent_line_nm": 253.639, "analyte_line_nm": 253.652,
         "separation_pm": 13, "severity": "medium",
         "note": "Co nearby line."},
    ],
    "se": [
        {"interferent": "Fe", "interferent_line_nm": 196.089, "analyte_line_nm": 195.995,
         "separation_pm": 94, "severity": "low",
         "note": "Fe nearby line at low wavelength region."},
        {"interferent": "Co", "interferent_line_nm": 196.026, "analyte_line_nm": 195.995,
         "separation_pm": 31, "severity": "low",
         "note": "Co nearby line."},
    ],
    "mn": [
        {"interferent": "Ga", "interferent_line_nm": 279.502, "analyte_line_nm": 279.482,
         "separation_pm": 20, "severity": "medium",
         "note": "Ga 279.502 close to Mn triplet line."},
    ],
    "ga": [
        {"interferent": "Fe", "interferent_line_nm": 287.417, "analyte_line_nm": 287.424,
         "separation_pm": 7, "severity": "high",
         "note": "Fe 287.417 very close to Ga 287.424."},
    ],
    "pb": [
        {"interferent": "Fe", "interferent_line_nm": 283.305, "analyte_line_nm": 283.306,
         "separation_pm": 1, "severity": "critical",
         "note": "Near exact overlap with Fe 283.305 — use Pb 217.0 nm as alternative."},
    ],
    "cr": [
        {"interferent": "Mn", "interferent_line_nm": 357.879, "analyte_line_nm": 357.869,
         "separation_pm": 10, "severity": "medium",
         "note": "Mn nearby line to Cr most sensitive line."},
        {"interferent": "V", "interferent_line_nm": 357.861, "analyte_line_nm": 357.869,
         "separation_pm": 8, "severity": "medium",
         "note": "V nearby line."},
    ],
    "ni": [
        {"interferent": "Fe", "interferent_line_nm": 232.037, "analyte_line_nm": 232.003,
         "separation_pm": 34, "severity": "low",
         "note": "Fe nearby line to Ni primary line."},
    ],
    "sb": [
        {"interferent": "Fe", "interferent_line_nm": 217.609, "analyte_line_nm": 217.581,
         "separation_pm": 28, "severity": "medium",
         "note": "Fe nearby line; use Sb 206.8 as alternative if needed."},
    ],
    "co": [
        {"interferent": "Hg", "interferent_line_nm": 240.727, "analyte_line_nm": 240.710,
         "separation_pm": 17, "severity": "low",
         "note": "Hg nearby line."},
    ],
    "as": [
        {"interferent": "Al", "interferent_line_nm": 193.759, "analyte_line_nm": 193.696,
         "separation_pm": 63, "severity": "low",
         "note": "Al molecular band near As line (hydride method recommended)."},
    ],
}


# ── 化学干扰数据库 ────────────────────────────────────────────
CHEMICAL_INTERFERENCE_DB: Dict[str, List[Dict[str, Any]]] = {
    # analyte -> list of (interferent, mechanism, effect, remedy)
    "ca": [
        {"interferent": "PO4^3-", "mechanism": "Formation of refractory Ca3(PO4)2",
         "effect": "Negative (depression)", "remedy": "Add LaCl3 or SrCl2 as releasing agent"},
        {"interferent": "Al", "mechanism": "Formation of refractory calcium aluminate",
         "effect": "Negative (depression)", "remedy": "Add LaCl3 or EDTA"},
        {"interferent": "SiO4^4-", "mechanism": "Formation of CaSiO3",
         "effect": "Negative (depression)", "remedy": "Add NH4Cl or HF digestion"},
        {"interferent": "SO4^2-", "mechanism": "Anion effect on atomization efficiency",
         "effect": "Variable", "remedy": "Matrix matching of standards"},
    ],
    "mg": [
        {"interferent": "Al", "mechanism": "Formation of MgAl2O4 spinel (refractory)",
         "effect": "Strong negative", "remedy": "Add Sr(NO3)2 or La(NO3)3 as releasing agent"},
        {"interferent": "Si", "mechanism": "Formation of refractory silicates",
         "effect": "Strong negative", "remedy": "Add Sr, La, or Ca as releasing agent"},
        {"interferent": "Ti", "mechanism": "Formation of mixed oxides",
         "effect": "Negative", "remedy": "Use N2O-C2H2 flame or add releasing agent"},
        {"interferent": "PO4^3-", "mechanism": "Phosphate binding",
         "effect": "Moderate negative", "remedy": "Add LaCl3"},
    ],
    "cr": [
        {"interferent": "Fe", "mechanism": "Intermetallic compound formation / competitive atomization",
         "effect": "Variable", "remedy": "Use rich air-C2H2 flame; add NH4Cl"},
        {"interferent": "Ni", "mechanism": "Competitive oxidation/reduction in flame",
         "effect": "Variable", "remedy": "Optimize fuel ratio"},
    ],
    "mo": [
        {"interferent": "Ca", "mechanism": "Formation of refractory CaMoO4",
         "effect": "Strong negative", "remedy": "Add Al, Sr, or use N2O-C2H2 rich flame"},
        {"interferent": "SO4^2-", "mechanism": "Sulfate suppression of Mo signal",
         "effect": "Negative", "remedy": "Convert to chloride form; add NH4NO3"},
    ],
    "ba": [
        {"interferent": "Al", "mechanism": "Formation of BaAl2O4",
         "effect": "Strong negative", "remedy": "N2O-C2H2 rich flame + releasing agent"},
        {"interferent": "PO4^3-", "mechanism": "Ba3(PO4)2 formation",
         "effect": "Strong negative", "remedy": "Add La or EDTA"},
    ],
    "sr": [
        {"interferent": "Al", "mechanism": "Refractory aluminate formation",
         "effect": "Strong negative", "remedy": "N2O-C2H2 + LaCl3"},
        {"interferent": "PO4^3-", "mechanism": "Sr3(PO4)2 formation",
         "effect": "Negative", "remedy": "Add La or Sr releasing agent"},
    ],
    "fe": [
        {"interferent": "Si", "mechanism": "Silicate formation reducing atomization",
         "effect": "Moderate negative", "remedy": "Add NH4Cl or CaCl2"},
    ],
    "na": [
        {"interferent": "PO4^3-", "mechanism": "Anion enhancement effect in cool flames",
         "effect": "Positive (enhancement)", "remedy": "Use leaner flame; matrix match standards"},
        {"interferent": "organic", "mechanism": "Organic solvent enhances nebulization efficiency",
         "effect": "Positive (enhancement)", "remedy": "Matrix matching or standard addition"},
    ],
    "k": [
        {"interferent": "organic", "mechanism": "Organic solvent enhances transport",
         "effect": "Positive (enhancement)", "remedy": "Matrix matching"},
    ],
}


# ── 易电离元素 (IEC) 数据库 ────────────────────────────────────
IONIZATION_EASY_ELEMENTS = {
    "li": {"ip_ev": 5.39, "suppressor": "KCl (1000 mg/L)"},
    "na": {"ip_ev": 5.14, "suppressor": "CsCl or KCl (1000 mg/L)"},
    "k":  {"ip_ev": 4.34, "suppressor": "CsCl (1000 mg/L)"},
    "rb": {"ip_ev": 4.18, "suppressor": "KCl or NaCl (1000 mg/L)"},
    "cs": {"ip_ev": 3.89, "suppressor": "NaCl or KCl (1000 mg/L)"},
    "ba": {"ip_ev": 5.21, "suppressor": "KCl (1000 mg/L)"},
    "sr": {"ip_ev": 5.69, "suppressor": "KCl (1000 mg/L)"},
    "ca": {"ip_ev": 6.11, "suppressor": "KCl (1000 mg/L)"},  # marginal in N2O-C2H2
}


@ChemMCPManager.register_tool
class AasInterferenceChecker(BaseTool):
    """
    AAS光谱干扰和化学干扰诊断工具。
    分析待测元素在给定基体中的潜在干扰，给出严重程度和消除建议。
    """
    __version__                = "0.1.0"
    name                       = "AasInterferenceChecker"
    func_name                  = "check_interferences"
    description                = ("Diagnose spectral, chemical, and ionization interferences "
                                 "in Atomic Absorption Spectroscopy (AAS) measurements.")
    implementation_description = ("Uses built-in databases of spectral line overlaps, chemical interference mechanisms, "
                                 "and ionization behavior to diagnose and recommend remedies for AAS interferences.")
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["AAS", "Interference", "Spectral Analysis",
                                   "Analytical Chemistry", "Quality Control"]
    required_envs              = []

    code_input_sig = [
        ("analyte_element",             "str",   "N/A",     "Analyte element symbol (e.g., 'Cd', 'Zn', 'Ca')."),
        ("matrix_elements",             "list",  "[]",      "List of potential interfering elements/ions in the sample matrix."),
        ("flame_type",                  "str",   "air-acetylene", "Flame type: 'air-acetylene' or 'nitrous_oxide-acetylene'."),
        ("wavelength_nm",               "float", "None",    "Analytical wavelength in nm (optional, for spectral overlap check)."),
        ("lamp_type",                   "str",   "HCL",     "Lamp type: 'HCL' (hollow cathode lamp) or 'EDL' (electrodeless discharge)."),
    ]

    text_input_sig = [
        ("input_params",                "str",   "N/A",     "Space-separated: analyte_element interferent1,interferent2,... [flame_type] [wavelength]"),
    ]

    output_sig = [
        ("spectral_interferences",      "list",  "List of spectral overlap interferences found."),
        ("chemical_interferences",       "list",  "List of chemical/matrix interferences found."),
        ("ionization_interference",      "dict",  "Ionization interference assessment (if applicable)."),
        ("recommended_remedies",         "list",  "Recommended actions to eliminate/minimize interferences."),
        ("severity_level",              "str",   "Overall severity: 'low', 'medium', 'high', or 'critical'."),
        ("summary",                     "str",   "Human-readable summary of the interference assessment."),
    ]

    examples = [
        {
            "code_input": {
                "analyte_element": "Cd",
                "matrix_elements": ["Fe"],
                "flame_type": "air-acetylene",
            },
            "text_input": {"input_params": "Cd Fe air-acetylene"},
            "output": {
                "severity_level": "critical",
                "summary": "",
            },
        },
        {
            "code_input": {
                "analyte_element": "Ca",
                "matrix_elements": ["PO4^3-", "Al"],
                "flame_type": "air-acetylene",
            },
            "text_input": {"input_params": "Ca PO4^3-,Al air-acetylene"},
            "output": {
                "severity_level": "high",
                "summary": "",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Load interference databases."""
        self.spectral_db = SPECTRAL_OVERLAP_DB
        self.chemical_db = CHEMICAL_INTERFERENCE_DB
        self.ionization_db = IONIZATION_EASY_ELEMENTS

    def _run_base(
        self,
        analyte_element: str,
        matrix_elements: Optional[List[str]] = None,
        flame_type: str = "air-acetylene",
        wavelength_nm: Optional[float] = None,
        lamp_type: str = "HCL",
    ) -> Dict[str, Any]:
        """Core diagnostic logic."""
        key = analyte_element.strip().lower()
        matrix = [m.strip().lower() for m in (matrix_elements or [])]

        spectral_intfs = []
        chemical_intfs = []
        remedies = set()
        severity_scores = []

        # ── 1. Spectral interference check ─────────────────────
        if key in self.spectral_db:
            for entry in self.spectral_db[key]:
                int_elem = entry["interferent"].lower()
                # Check if this interferent is in the provided matrix
                is_in_matrix = not matrix or int_elem in matrix or any(
                    int_elem in m or m in int_elem for m in matrix
                )
                if is_in_matrix:
                    spectral_intfs.append(entry)
                    sev = entry.get("severity", "unknown")
                    severity_scores.append({"type": "spectral", "severity": sev})
                    note = entry.get("note", "")
                    if sev == "critical":
                        remedies.add(f"CRITICAL: {note} Use alternative analytical wavelength or high-resolution mode.")
                    elif sev == "high":
                        remedies.add(f"HIGH: {note} Consider Zeeman background correction or alternate line.")

        # ── 2. Chemical interference check ─────────────────────
        if key in self.chemical_db:
            for entry in self.chemical_db[key]:
                int_chem = entry["interferent"].lower()
                is_in_matrix = not matrix or int_chem in matrix or any(
                    int_chem in m or m.replace("^", "").replace("+", "").replace("-", "") in int_chem.replace("^", "").replace("+", "").replace("-", "")
                    for m in matrix
                )
                if is_in_matrix:
                    chemical_intfs.append(entry)
                    remedies.add(entry.get("remedy", ""))
                    severity_scores.append({"type": "chemical", "severity": "high"})

        # ── 3. Ionization interference check ───────────────────
        ion_info = None
        if key in self.ionization_db:
            entry = self.ionization_db[key]
            ip = entry["ip_ev"]
            # In N2O-C2H2 flame (~2950 K), elements with IP < ~7 eV may ionize
            # In air-C2H2 flame (~2300 K), elements with IP < ~6 eV may ionize
            is_n2o = "nitrous" in flame_type.lower()
            threshold_ip = 7.0 if is_n2o else 6.0

            if ip < threshold_ip:
                ion_info = {
                    "element": key.upper(),
                    "ionization_potential_eV": ip,
                    "at_risk": True,
                    "flame_type": flame_type,
                    "recommended_suppressor": entry["suppressor"],
                    "explanation": (
                        f"{key.upper()} has low ionization potential ({ip} eV). "
                        f"In {flame_type} flame, significant ionization occurs, "
                        f"reducing ground-state atom population and sensitivity."
                    ),
                }
                remedies.add(f"Add ionization suppressor: {entry['suppressor']}")
                severity_scores.append({"type": "ionization", "severity": "medium"})
            else:
                ion_info = {
                    "element": key.upper(),
                    "ionization_potential_eV": ip,
                    "at_risk": False,
                    "flame_type": flame_type,
                    "recommended_suppressor": entry["suppressor"],
                    "explanation": f"Ionization risk is low in {flame_type} flame.",
                }

        # ── 4. Determine overall severity ──────────────────────
        severity_order = {"low": 1, "medium": 2, "high": 3, "critical": 4}
        max_sev = "low"
        for s in severity_scores:
            sv = s.get("severity", "low")
            if severity_order.get(sv, 0) > severity_order.get(max_sev, 0):
                max_sev = sv

        # Add general remedies
        if spectral_intfs:
            remedies.add("Consider using Zeeman or D₂ background correction for spectral overlaps.")
        if chemical_intfs:
            remedies.add("Standard addition method can compensate for many chemical interferences.")
        if not spectral_intfs and not chemical_intfs and (not ion_info or not ion_info["at_risk"]):
            remedies.add("No significant interferences detected under specified conditions.")

        # Build summary
        parts = [f"AAS interference analysis for **{key.upper()}** ({flame_type}):"]
        if spectral_intfs:
            parts.append(f"  - Spectral interferences: {len(spectral_intfs)} found")
        if chemical_intfs:
            parts.append(f"  - Chemical interferences: {len(chemical_intfs)} found")
        if ion_info and ion_info["at_risk"]:
            parts.append(f"  - Ionization interference: YES (IP={ion_info['ionization_potential_eV']} eV)")
        elif ion_info:
            parts.append("  - Ionization interference: not significant")
        parts.append(f"  - Overall severity: **{max_sev.upper()}**")

        logger.info(f"AAS interference check for {key.upper()}: severity={max_sev}, "
                     f"spectral={len(spectral_intfs)}, chemical={len(chemical_intfs)}")
        return {
            "spectral_interferences": spectral_intfs,
            "chemical_interferences": chemical_intfs,
            "ionization_interference": ion_info,
            "recommended_remedies": sorted(remedies),
            "severity_level": max_sev,
            "summary": "\n".join(parts),
        }

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """Parse text input."""
        try:
            parts = input_params.split()
            if not parts:
                raise ValueError("Empty input.")

            analyte = parts[0]
            matrix = None
            flame = "air-acetylene"
            wavelength = None

            idx = 1
            while idx < len(parts):
                p = parts[idx]
                if p.endswith(",") or any(c.isalpha() for c in p):
                    if matrix is None:
                        matrix = [x.strip() for x in p.split(",") if x.strip()]
                    else:
                        matrix.extend([x.strip() for x in p.split(",") if x.strip()])
                else:
                    try:
                        float(p)
                        if wavelength is None:
                            wavelength = float(p)
                    except ValueError:
                        flame = p
                idx += 1

            return self._run_base(analyte, matrix, flame, wavelength)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
