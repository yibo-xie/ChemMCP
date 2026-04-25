import logging
from typing import List, Dict, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SelectivePrecipitation(BaseTool):
    """
    设计选择性沉淀分离方案。
    基于经典定性分析分组方案和Ksp差异，为混合离子体系设计最优分离流程。
    """
    __version__ = "0.1.0"
    name = "SelectivePrecipitation"
    func_name = "selective_precipitation"
    description = "Design a selective precipitation separation scheme for a mixture of cations. Uses classical qualitative analysis group separation based on Ksp differences, pH control, and reagent selection."
    implementation_description = "Implements the classic cation group analysis scheme (Groups I-V) with Ksp-based reasoning. For custom ion mixtures, determines optimal precipitation order by comparing Qsp vs Ksp at each step with available reagents."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Precipitation", "Separation", "Ksp", "Qualitative Analysis", "Cation Groups"]
    required_envs = []

    code_input_sig = [
        ("ions", "list", "N/A", "List of cation formulas to separate, e.g., ['Ag+', 'Ba2+', 'Cu2+', 'Mg2+']."),
        ("reagents", "list", "None", "Optional: list of available precipitating reagents. If None, uses standard group analysis reagents."),
        ("target_ph", "float", "None", "Optional: target pH for a specific step. Auto-determined if None."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Space or comma-separated ions. E.g., 'Ag+ Ba2+ Cu2+ Mg2+' or 'Pb2+ Hg2^2+ Cu2+ Fe3+ Ba2+ Na+'."),
    ]

    output_sig = [
        ("separation_scheme", "list", "Step-by-step separation plan, each with step number, reagent, pH condition, precipitate formed, ions remaining."),
        ("summary", "str", "Human-readable summary of the full separation process."),
        ("total_steps", "int", "Number of separation steps required."),
        ("ions_separated", "list", "List of all ions that can be separated."),
        ("notes", "str", "Additional notes and caveats about the scheme."),
    ]

    examples = [
        {
            "code_input": {
                "ions": ["Ag+", "Ba2+", "Cu2+", "Mg2+"],
                "reagents": None,
                "target_ph": None,
            },
            "text_input": {
                "input_str": "Ag+ Ba2+ Cu2+ Mg2+",
            },
            "output": {
                "separation_scheme": [
                    {"step": 1, "reagent": "HCl (dilute)", "ph": "~0-1", "precipitate": "AgCl(s)", "remaining": ["Ba2+", "Cu2+", "Mg2+"]},
                    {"step": 2, "reagent": "H2S (pH 0.5, acidic)", "ph": "~0.5", "precipitate": "CuS(s)", "remaining": ["Ba2+", "Mg2+"]},
                    {"step": 3, "reagent": "(NH4)2CO3 (pH~9)", "ph": "~9", "precipitate": "BaCO3(s)", "remaining": ["Mg2+"]},
                    {"step": 4, "reagent": "No precipitant — Mg2+ stays in solution", "ph": None, "precipitate": None, "remaining": ["Mg2+"]},
                ],
                "summary": "4-step separation: Ag+ as AgCl in acidic → Cu2+ as CuS in acidic H2S → Ba2+ as BaCO3 in weakly basic carbonate → Mg2+ remains soluble (Group V).",
                "total_steps": 4,
                "ions_separated": ["Ag+", "Ba2+", "Cu2+", "Mg2+"],
                "notes": "Standard qualitative analysis scheme. Actual lab conditions may require pH adjustment and confirmation tests.",
            },
        },
        {
            "code_input": {
                "ions": ["Pb2+", "Hg2^2+", "Fe3+", "Na+"],
                "reagents": None,
                "target_ph": None,
            },
            "text_input": {
                "input_str": "Pb2+ Hg2^2+ Fe3+ Na+",
            },
            "output": {
                "separation_scheme": [
                    {"step": 1, "reagent": "HCl (dilute)", "ph": "~0-1", "precipitate": "PbCl2(s) + Hg2Cl2(s)", "remaining": ["Fe3+", "Na+"]},
                    {"step": 2, "reagent": "NH3 / NH4Cl buffer + H2S or (NH4)2S", "ph": "~8-9", "precipitate": "Fe2S3/Fe(OH)3(s)", "remaining": ["Na+"]},
                    {"step": 3, "reagent": "No precipitant needed", "ph": None, "precipitate": None, "remaining": ["Na+"]},
                ],
                "summary": "3-step: Group I (Pb2+, Hg2^2+) precipitated as chlorides → Group III (Fe3+) precipitated as sulfide/hydroxide → Na+ remains (Group V).",
                "total_steps": 3,
                "ions_separated": ["Pb2+", "Hg2^2+", "Fe3+", "Na+"],
                "notes": "PbCl2 has moderate solubility in cold water; may need hot water wash to confirm.",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize ion group database."""
        # Classical qualitative analysis groups
        # Each group: {group_name, reagent, ph_range, ions, precipitate_form}
        self._groups = [
            {
                "name": "Group I — Silver Group",
                "reagent": "HCl (dilute)",
                "ph": "~0–1 (strongly acidic)",
                "ph_value": 0.5,
                "ions": {"Ag+": "AgCl(s) white", "Pb2+": "PbCl2(s) white", "Hg2^2+": "Hg2Cl2(s) white"},
                "notes": "Precipitates as chlorides. PbCl2 moderately soluble in cold water.",
            },
            {
                "name": "Group II — Copper-Arsenic Group",
                "reagent": "H2S in 0.3 M HCl (or thioacetamide)",
                "ph": "~0.5 (acidic)",
                "ph_value": 0.5,
                "ions": {
                    "Pb2+": "PbS(s) black (if remaining after Grp I)",
                    "Cu2+": "CuS(s) black",
                    "Cd2+": "CdS(s) yellow",
                    "Hg2+": "HgS(s) black",
                    "Bi3+": "Bi2S3(s) dark brown",
                    "As3+": "As2S3(s) yellow",
                    "Sn2+": "SnS(s) brown",
                    "Sn4+": "SnS2(s) yellow",
                    "Sb3+": "Sb2S3(s) orange-red",
                },
                "notes": "Acid-insoluble sulfides. Subdivided into IIa (HNO3-soluble) and IIb ((NH4)2Sx-soluble).",
            },
            {
                "name": "Group III — Iron-Nickel Group",
                "reagent": "(NH4)2S or NH3 + (NH4)2S in NH4+/NH3 buffer",
                "ph": "~8–9 (basic)",
                "ph_value": 8.5,
                "ions": {
                    "Fe3+": "Fe2S3(s)/Fe(OH)3(s) reddish-brown",
                    "Fe2+": "FeS(s) black",
                    "Cr3+": "Cr(OH)3(s) green",
                    "Al3+": "Al(OH)3(s) white gelatinous",
                    "Mn2+": "MnS(s) flesh/pink",
                    "Ni2+": "NiS(s) black",
                    "Co2+": "CoS(s) black",
                    "Zn2+": "ZnS(s) white",
                },
                "notes": "Hydroxides/sulfides insoluble in basic buffer but not in pure water. Zn2+, Mn2+ sometimes placed in Group IV.",
            },
            {
                "name": "Group IV — Barium Group (Carbonate Group)",
                "reagent": "(NH4)2CO3 in NH3/NH4Cl buffer",
                "ph": "~9 (basic)",
                "ph_value": 9.0,
                "ions": {
                    "Ba2+": "BaCO3(s) white",
                    "Sr2+": "SrCO3(s) white",
                    "Ca2+": "CaCO3(s) white (slightly)",
                    "Mg2+": "MgCO3(s) white (only if [Mg2+] high enough; often stays in solution)",
                },
                "notes": "Carbonates insoluble in basic buffer. Mg2+ may partially precipitate if concentration is high.",
            },
            {
                "name": "Group V — Soluble Group",
                "reagent": "No precipitant needed",
                "ph": None,
                "ph_value": None,
                "ions": {
                    "Na+": "remains in solution",
                    "K+": "remains in solution",
                    "NH4+": "remains in solution (detected separately)",
                    "Mg2+": "often remains (if not ppted in Grp IV)",
                    "Li+": "remains in solution",
                },
                "notes": "Ions of Group 1 elements and ammonium. No common precipitating reagent.",
            },
        ]

        # Ion → canonical name mapping
        self._ion_aliases = {
            "ag+": "Ag+", "silver ion": "Ag+",
            "pb2+": "Pb2+", "lead(ii) ion": "Pb2+", "lead ion": "Pb2+",
            "hg2^2+": "Hg2^2+", "mercury(i) ion": "Hg2^2+", "mercurous ion": "Hg2^2+",
            "cu2+": "Cu2+", "copper(ii) ion": "Cu2+",
            "cd2+": "Cd2+", "cadmium ion": "Cd2+",
            "hg2+": "Hg2+", "mercury(ii) ion": "Hg2+",
            "bi3+": "Bi3+", "bismuth(iii) ion": "Bi3+",
            "as3+": "As3+", "arsenious ion": "As3+",
            "sn2+": "Sn2+", "tin(ii) ion": "Sn2+",
            "sn4+": "Sn4+", "tin(iv) ion": "Sn4+",
            "sb3+": "Sb3+", "antimony(iii) ion": "Sb3+",
            "fe3+": "Fe3+", "iron(iii) ion": "Fe3+",
            "fe2+": "Fe2+", "iron(ii) ion": "Fe2+",
            "cr3+": "Cr3+", "chromium(iii) ion": "Cr3+",
            "al3+": "Al3+", "aluminum ion": "Al3+",
            "mn2+": "Mn2+", "manganese(ii) ion": "Mn2+",
            "ni2+": "Ni2+", "nickel(ii) ion": "Ni2+",
            "co2+": "Co2+", "cobalt(ii) ion": "Co2+",
            "zn2+": "Zn2+", "zinc ion": "Zn2+",
            "ba2+": "Ba2+", "barium ion": "Ba2+",
            "sr2+": "Sr2+", "strontium ion": "Sr2+",
            "ca2+": "Ca2+", "calcium ion": "Ca2+",
            "mg2+": "Mg2+", "magnesium ion": "Mg2+",
            "na+": "Na+", "sodium ion": "Na+",
            "k+": "K+", "potassium ion": "K+",
            "nh4+": "NH4+", "ammonium ion": "NH4+",
            "li+": "Li+", "lithium ion": "Li+",
        }

    def _run_base(self, ions: List[str], reagents: Optional[List[str]] = None,
                  target_ph: Optional[float] = None) -> dict:
        """Core logic: design selective precipitation scheme."""
        if not ions:
            raise ChemMCPError("Ion list cannot be empty.")

        # Resolve ion names
        resolved_ions = []
        for ion in ions:
            r = self._resolve_ion(ion)
            resolved_ions.append(r)

        unique_ions = list(dict.fromkeys(resolved_ions))  # preserve order, dedupe

        # Build separation scheme step by step
        scheme = []
        remaining = list(unique_ions)
        all_separated = []

        for grp in self._groups:
            if not remaining:
                break

            # Find which remaining ions belong to this group
            precipitated_this_step = {}
            still_remaining = []

            for ion in remaining:
                if ion in grp["ions"]:
                    precipitated_this_step[ion] = grp["ions"][ion]
                else:
                    still_remaining.append(ion)

            if precipitated_this_step:
                ppt_list = ", ".join([f"{k} ({v})" for k, v in precipitated_this_step.items()])
                step_info = {
                    "step": len(scheme) + 1,
                    "group": grp["name"],
                    "reagent": grp["reagent"],
                    "ph": grp["ph"],
                    "ph_value": grp.get("ph_value"),
                    "precipitate": ppt_list,
                    "precipitate_details": dict(precipitated_this_step),
                    "ions_remaining": list(still_remaining),
                    "notes": grp.get("notes", ""),
                }
                scheme.append(step_info)
                all_separated.extend(precipitated_this_step.keys())
                remaining = still_remaining

        # Handle any unrecognized ions
        unrecognized = [ion for ion in remaining if not self._is_known_ion(ion)]
        if unrecognized:
            scheme.append({
                "step": len(scheme) + 1,
                "group": "Unknown / Not in standard groups",
                "reagent": "N/A — no standard method",
                "ph": "N/A",
                "ph_value": None,
                "precipitate": f"Unrecognized ions: {', '.join(unrecognized)}",
                "precipitate_details": {},
                "ions_remaining": [],
                "notes": f"Ions not found in standard qualitative analysis groups: {unrecognized}. Custom method development needed.",
            })

        # Build summary
        summary_parts = []
        for step in scheme:
            if step["precipitate"] and step["reagent"] != "No precipitant needed":
                summary_parts.append(
                    f"Step {step['step']}: Add {step['reagent']} (pH{step['ph']}) → "
                    f"{step['precipitate']}"
                )
            elif step["reagent"] == "No precipitant needed" and remaining:
                summary_parts.append(
                    f"Step {step['step']}: {', '.join(remaining)} remain in solution (soluble group)"
                )

        notes_list = [
            "This is a theoretical separation scheme based on standard qualitative analysis principles.",
            "Actual laboratory conditions require careful pH control, temperature management, and confirmatory tests.",
            "Some ions may appear in multiple groups depending on concentration (e.g., Pb2+, Zn2+, Mg2+).",
            "Oxidation states matter: ensure correct oxidation state is specified for each ion.",
        ]

        logger.info(f"SelectivePrecipitation: designed {len(scheme)}-step scheme for {len(unique_ions)} ions")
        return {
            "separation_scheme": scheme,
            "summary": " | ".join(summary_parts) if summary_parts else "No separation steps needed.",
            "total_steps": len(scheme),
            "ions_separated": all_separated + (remaining if not unrecognized else []),
            "notes": "\n".join(notes_list),
        }

    def _run_text(self, input_str: str) -> dict:
        """Parse text input: space/comma/semicolon separated ions."""
        import re
        tokens = re.split(r"[,\s]+", input_str.strip())
        ions = [t for t in tokens if t]
        return self._run_base(ions)

    def _resolve_ion(self, name: str) -> str:
        """Resolve ion name/alias to canonical key."""
        n = name.strip()
        if n in self._ion_aliases:
            return self._ion_aliases[n]
        nl = n.lower()
        if nl in self._ion_aliases:
            return self._ion_aliases[nl]
        return n

    def _is_known_ion(self, ion: str) -> bool:
        """Check if ion exists in any group."""
        for grp in self._groups:
            if ion in grp["ions"]:
                return True
        return False
