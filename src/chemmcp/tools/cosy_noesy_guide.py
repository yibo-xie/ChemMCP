"""
COSY/NOESY Guide - 2D NMR interpretation guide for through-bond (COSY)
and through-space (NOESY) connectivity analysis.
"""

import logging
from typing import Dict, List, Tuple, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CosyNoesyGuide(BaseTool):
    __version__      = "0.1.0"
    name             = "CosyNoesyGuide"
    func_name        = "guide_2d_nmr"
    description      = "2D NMR (COSY/NOESY) interpretation guide: analyze cross-peaks to map spin systems and spatial proximity."
    implementation_description = "Parses 1D NMR peaks and 2D cross-peak correlations to identify spin systems (COSY: J-coupled neighbors) or spatial relationships (NOESY: through-space <5Å)."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Molecule"]
    tags             = ["2D NMR", "COSY", "NOESY", "Spectroscopy", "Spin System"]
    required_envs    = []

    code_input_sig   = [
        ("peaks_1d", "list", "N/A", "List of 1D ¹H NMR peaks as [(chemical_shift_ppm, multiplicity, integration), ...]."),
        ("cross_peaks", "list", "[]", "List of 2D cross-peaks as [(f1_shift, f2_shift), ...] or [(f1_shift, f2_shift, intensity), ...]."),
        ("experiment_type", "str", "cosy", "Experiment type: 'cosy' for COSY (through-bond coupling) or 'noesy' for NOESY (through-space)."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Pipe-separated: '1D_peaks|cross_peaks|experiment_type'. Example: '(1.0,t,3H),(1.4,m,2H)|(1.0,1.4)|cosy'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict containing: experiment_type, spin_systems (list of connected proton groups), connectivity_map, structural_inferences, recommended_next_experiments."),
    ]

    examples         = [
        {
            "code_input": {
                "peaks_1d": [[7.25, "m", 5], [4.15, "q", 2], [3.65, "s", 3], [2.85, "t", 2], [1.35, "t", 3]],
                "cross_peaks": [[4.15, 2.85], [2.85, 1.35]],
                "experiment_type": "cosy",
            },
            "text_input": {"input_params": "(7.25,m,5),(4.15,q,2),(3.65,s,3),(2.85,t,2),(1.35,t,3)|(4.15,2.85),(2.85,1.35)|cosy"},
            "output": {
                "result": {
                    "experiment_type": "COSY",
                    "spin_systems": [
                        {"members": [4.15, 2.85, 1.35], "type": "coupled chain"},
                        {"members": [7.25], "type": "isolated"},
                        {"members": [3.65], "type": "isolated"},
                    ],
                    "connectivity_map": {
                        "4.15 ↔ 2.85": "J-coupled (COSY correlation)",
                        "2.85 ↔ 1.35": "J-coupled (COSY correlation)",
                    },
                    "structural_inferences": ["Ethyl-like fragment: CH(q,2H)-CH₂(t,2H)-CH₃(t,3H)", "Aromatic protons at 7.25 ppm (5H)", "Isolated singlet at 3.65 ppm (3H) - likely OCH₃"],
                    "recommended_next_experiments": ["HSQC to correlate H-C pairs", "HMBC for long-range C-H connectivities", "NOESY if stereochemistry is of interest"],
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, peaks_1d: list, cross_peaks: list = None, experiment_type: str = "cosy") -> dict:
        """Core logic: analyze 2D NMR correlations."""
        if not peaks_1d:
            raise ChemMCPError("1D NMR peak list is required.")

        exp_type = experiment_type.lower().strip()
        if exp_type not in ("cosy", "noesy"):
            raise ChemMCPError("experiment_type must be 'cosy' or 'noesy'")

        cross_peaks = cross_peaks or []

        # Parse 1D peaks into structured format
        parsed_1d = []
        for p in peaks_1d:
            shift = round(float(p[0]), 2)
            mult = str(p[1]) if len(p) > 1 else "unknown"
            integ = float(p[2]) if len(p) > 2 else 0
            parsed_1d.append({"shift": shift, "multiplicity": mult, "integration": integ})

        all_shifts = set(pp["shift"] for pp in parsed_1d)

        # Parse cross-peaks
        correlations: List[Dict[str, Any]] = []
        for cp in cross_peaks:
            f1 = round(float(cp[0]), 2)
            f2 = round(float(cp[1]), 2)
            intensity = str(cp[2]) if len(cp) > 2 else "medium"

            if f1 not in all_shifts and f2 not in all_shifts:
                continue  # skip if neither dimension matches a known peak

            corr_type = "J-coupling (³JHH)" if exp_type == "cosy" else "Through-space dipolar (<5Å)"
            correlations.append({
                "f1": f1,
                "f2": f2,
                "intensity": intensity,
                "interpretation": f"{f1} ↔ {f2}: {corr_type}",
            })

        # Build spin systems using union-find on correlations
        spin_systems = self._build_spin_systems(all_shifts, correlations)

        # Generate structural inferences
        inferences = self._infer_structure(parsed_1d, spin_systems, correlations, exp_type)

        # Recommend next experiments
        next_exp = []
        if exp_type == "cosy":
            next_exp.extend(["HSQC (¹H-¹³C one-bond correlation)", "HMBC (²J/³J CH long-range)", "NOESY / ROESY for spatial information"])
        else:
            next_exp.extend(["COSY for through-bond connectivity", "ROESY for medium-sized molecules"])

        return {
            "experiment_type": exp_type.upper(),
            "spin_systems": spin_systems,
            "connectivity_map": {c["interpretation"].split(":")[0]: c["interpretation"].split(": ", 1)[1] for c in correlations},
            "structural_inferences": inferences,
            "recommended_next_experiments": next_exp,
        }

    @staticmethod
    def _build_spin_systems(all_shifts: set, correlations: list) -> list:
        """Group shifts into spin systems via union-find on correlations."""
        parent = {s: s for s in all_shifts}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[rx] = ry

        for c in correlations:
            if c["f1"] in parent and c["f2"] in parent:
                union(c["f1"], c["f2"])

        systems_dict: Dict[str, List[float]] = {}
        for s in all_shifts:
            root = find(s)
            if root not in systems_dict:
                systems_dict[root] = []
            systems_dict[root].append(s)

        result = []
        for members in systems_dict.values():
            stype = "coupled chain" if len(members) > 1 else "isolated"
            result.append({
                "members": sorted(members),
                "size": len(members),
                "type": stype,
            })

        return sorted(result, key=lambda x: (-x["size"], x["members"][0]))

    @staticmethod
    def _infer_structure(parsed_1d: list, spin_systems: list, correlations: list, exp_type: str) -> list:
        """Generate structural interpretations from the data."""
        inferences = []

        for ss in spin_systems:
            members = ss["members"]
            member_data = [(p["shift"], p["multiplicity"], p["integration"]) for p in parsed_1d if p["shift"] in members]

            if len(members) >= 3:
                # Look for ethyl pattern: quartet + triplet
                q_signals = [m for m in member_data if m[1] in ("q", "quartet")]
                t_signals = [m for m in member_data if m[1] in ("t", "triplet")]
                if q_signals and t_signals:
                    inferences.append(f"Ethyl group detected: CH at {q_signals[0][0]} ppm ({q_signals[0][2]}H) coupled to CH₂ at {t_signals[0][0]} ppm ({t_signals[0][2]}H)")

                # Isopropyl pattern: doublet + septet/heptet
                d_signals = [m for m in member_data if m[1] in ("d", "doublet")]
                septet_signals = [m for m in member_data if "septet" in m[1] or "heptet" in m[1]]
                if d_signals and len(d_signals) >= 2 and septet_signals:
                    total_d_h = sum(d[2] for d in d_signals)
                    inferences.append(f"Isopropyl group likely: {len(d_signals)} doublets ({total_d_h}H total) + septet at {septet_signals[0][0]} ppm ({septet_signals[0][2]}H)")

            elif len(members) == 2:
                inferences.append(f"Two-spin system: δ {member_data[0][0]} ({member_data[0][1]}, {member_data[0][2]}H) ↔ δ {member_data[1][0]} ({member_data[1][1]}, {member_data[1][2]}H)")

        # Check for characteristic isolated signals
        for p in parsed_1d:
            shift = p["shift"]
            mult = p["multiplicity"]
            integ = p["integration"]

            # Check if this shift is isolated (in a size-1 spin system)
            is_isolated = any(ss["size"] == 1 and shift in ss["members"] for ss in spin_systems)

            if is_isolated:
                if mult == "s":
                    if 6.5 <= shift <= 8.5:
                        inferences.append(f"Aromatic/isolated olefinic singlet at {shift} ppm ({integ}H)")
                    elif 3.3 <= shift <= 4.5:
                        inferences.append(f"Isolated singlet at {shift} ppm ({integ}H): possibly exchangeable (OH/NH) or OCH₃/Si(CH₃)₃")
                    elif 9.0 <= shift <= 10.5:
                        inferences.append(f"Aldehyde proton at {shift} ppm ({integ}H)")
                    elif shift < 2.0:
                        inferences.append(f"Aliphatic singlet at {shift} ppm ({integ}H): possibly tert-butyl or equivalent group")
                elif mult == "d" and integ == 6:
                    inferences.append(f"Doublet at {shift} ppm ({integ}H): two overlapping methyl groups (e.g., isopropyl)")

        # NOESY-specific inferences
        if exp_type == "noesy" and correlations:
            for c in correlations:
                f1, f2 = c["f1"], c["f2"]
                inferences.append(f"NOE between {f1} and {f2} ppm → these protons are within ~5 Å in space")

        return inferences

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split("|")
            peaks_str = parts[0].strip() if len(parts) > 0 else ""
            cross_str = parts[1].strip() if len(parts) > 1 else ""
            exp_type = parts[2].strip().lower() if len(parts) > 2 else "cosy"

            # Parse 1D peaks
            peaks_1d = []
            for m in __import__("re").findall(r'\(([^)]+)\)', peaks_str):
                items = [x.strip() for x in m.split(",")]
                if len(items) >= 3:
                    peaks_1d.append([float(items[0]), items[1], float(items[2])])

            # Parse cross-peaks
            cross_peaks = []
            for m in __import__("re").findall(r'\(([^)]+)\)', cross_str):
                items = [x.strip() for x in m.split(",")]
                if len(items) >= 2:
                    cross_peaks.append([float(items[0]), float(items[1])] + ([items[2]] if len(items) > 2 else []))

            return self._run_base(peaks_1d, cross_peaks, exp_type)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: '1D_peaks|cross_peaks|experiment_type'")
