"""
光谱化学序列查询工具
Spectrochemical series query and comparison for coordination chemistry.
"""
import logging
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SpectrochemicalSeries(BaseTool):
    """
    查询和比较光谱化学序列。
    提供完整配体场强序列、金属离子趋势、Δo/Δt比值等信息。
    """
    __version__ = "0.1.0"
    name = "SpectrochemicalSeries"
    func_name = "spectrochemical_series"
    description = "Query and compare spectrochemical series: ligand field strength ordering, metal ion trends, Δo/Δt ratios, and splitting parameter estimates."
    implementation_description = "Provides complete spectrochemical series data including ligand series (weak → strong), metal oxidation state trends (2+ < 3+ < 4+), period trends (3d < 4d < 5d), and relative Δo values with literature references."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Coordination Chemistry", "Spectrochemical Series", "Crystal Field", "Ligand Field", "Splitting"]
    required_envs = []

    code_input_sig = [
        ("query_type", "str", "full_series", "Type of query: 'full_series', 'compare_ligands', 'metal_trend', 'ligand_position', 'delta_estimate'."),
        ("ligands", "list", "None", "List of ligand names to compare or look up (optional)."),
        ("metal_ion", "str", "None", "Metal ion for metal-specific queries (optional)."),
        ("geometry", "str", "octahedral", "'octahedral' or 'tetrahedral'."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query string: 'query_type [ligand1 ligand2 ...] [metal_ion]', e.g., 'compare_ligands NH3 H2O CN-' or 'metal_trend Fe'."),
    ]

    output_sig = [
        ("query_type", "str", "The type of query performed."),
        ("series_data", "dict", "The spectrochemical series data requested."),
        ("explanation", "str", "Explanation of the data and its chemical significance."),
        ("references", "list", "Literature references for the data."),
    ]

    examples = [
        {
            "code_input": {
                "query_type": "full_series",
                "ligands": None,
                "metal_ion": None,
                "geometry": "octahedral",
            },
            "text_input": {
                "query": "full_series"
            },
            "output": {
                "query_type": "full_series",
                "series_data": {
                    "ligand_series": "I⁻ < Br⁻ < S²⁻ < SCN⁻ < Cl⁻ < NO₃⁻ < N₃⁻ < F⁻ < OH⁻ < oxalate²⁻ < H₂O < NCS⁻ < CH₃CN < py < NH₃ < en < bipy < phen < NO₂⁻ < PPh₃ < CN⁻ < CO",
                    "metal_trend": "Mn²⁺ < Ni²⁺ < Co²⁺ < Fe²⁺ < V²⁺ < Fe³⁺ < Co³⁺ < Cr³⁺ < Os³⁺ < Ir³⁺ < Pt⁴⁺ < Pd³⁺",
                    "oxidation_state": "higher oxidation state → larger Δo",
                    "period_trend": "3d (first row) < 4d (second row) < 5d (third row)",
                    "delta_t_ratio": "Δt ≈ (4/9) × Δo for same metal/ligands",
                },
                "explanation": "Full spectrochemical series showing ligand field strength from weakest to strongest...",
                "references": ["Figgis & Hitchman (2000)", "Huheey et al. (1993)", "Shriver & Atkins (2010)"],
            }
        },
        {
            "code_input": {
                "query_type": "compare_ligands",
                "ligands": ["NH3", "H2O", "CN-", "Cl-"],
                "metal_ion": None,
                "geometry": "octahedral",
            },
            "text_input": {
                "query": "compare_ligands NH3 H2O CN- Cl-"
            },
            "output": {
                "query_type": "compare_ligands",
                "series_data": {
                    "ordering": "Cl⁻ (weakest) < H₂O < NH₃ < CN⁻ (strongest)",
                    "relative_delta": {"Cl-": 0.80, "H2O": 1.00, "NH3": 1.25, "CN-": 1.70},
                    "field_category": {"Cl-": "weak", "H2O": "weak/intermediate", "NH3": "intermediate/strong", "CN-": "strong"},
                },
                "explanation": "Comparison of ligand field strengths with respect to H2O as reference (1.0)...",
                "references": [],
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize spectrochemical series database."""
        # Full ligand spectrochemical series (weak field → strong field)
        # Each entry: (name, symbol, relative_factor vs H2O=1.0, category)
        self._ligand_series = [
            ("iodide", "I⁻", 0.70, "very weak"),
            ("bromide", "Br⁻", 0.76, "very weak"),
            ("sulfide", "S²⁻", 0.80, "very weak"),
            ("thiocyanato-S", "SCN⁻ (S)", 0.76, "weak"),
            ("chloride", "Cl⁻", 0.80, "weak"),
            ("nitrate", "NO₃⁻", 0.85, "weak"),
            ("azide", "N₃⁻", 0.87, "weak"),
            ("fluoride", "F⁻", 0.90, "weak"),
            ("hydroxide", "OH⁻", 0.83, "weak"),
            ("oxalate", "C₂O₄²⁻", 0.99, "weak-intermediate"),
            ("water", "H₂O", 1.00, "reference (intermediate)"),
            ("thiocyanato-N", "NCS⁻ (N)", 1.05, "intermediate"),
            ("acetonitrile", "CH₃CN", 1.10, "intermediate"),
            ("pyridine", "py", 1.18, "intermediate-strong"),
            ("ammonia", "NH₃", 1.25, "intermediate-strong"),
            ("ethylenediamine", "en", 1.28, "strong"),
            ("bipyridine", "bipy", 1.30, "strong"),
            ("phenanthroline", "phen", 1.34, "strong"),
            ("nitro", "NO₂⁻", 1.50, "strong"),
            ("triphenylphosphine", "PPh₃", 1.65, "strong"),
            ("cyanide", "CN⁻", 1.70, "very strong"),
            ("carbonyl", "CO", 1.85, "very strong"),
        ]

        # Metal ion trend (increasing Δo for same ligand)
        self._metal_trend = [
            ("Mn²⁺", 0.75, "3d⁵, weak"), ("Ni²⁺", 0.82, "3d⁸"),
            ("Co²⁺", 0.89, "3d⁷"), ("Fe²⁺", 0.95, "3d⁶"),
            ("V²⁺", 1.05, "3d³"), ("Fe³⁺", 1.15, "3d⁵"),
            ("Co³⁺", 1.35, "3d⁶ low-spin"), ("Cr³⁺", 1.45, "3d³"),
            ("Os³⁺", 1.65, "5d⁵"), ("Ir³⁺", 1.75, "5d⁶"),
            ("Pt⁴⁺", 1.90, "5d⁶"), ("Pd³⁺", 1.80, "4d⁸"),
            ("Rh³⁺", 1.55, "4d⁶"), ("Ru³⁺", 1.50, "4d⁵"),
            ("Mo³⁺", 1.30, "4d³"),
        ]

        # Period comparison: Δo(4d) ≈ 1.3–1.5 × Δo(3d); Δo(5d) ≈ 1.5–2.0 × Δo(3d)
        self._period_factors = {
            "3d": 1.0, "4d": 1.4, "5d": 1.8,
        }

        # Oxidation state effect: Δo increases ~50-100% per unit increase in OS
        self._ox_factors = {"+1": 0.6, "+2": 1.0, "+3": 1.5, "+4": 2.0}

    def _run_base(self, query_type: str = "full_series", ligands: list = None,
                  metal_ion: str = None, geometry: str = "octahedral") -> dict:
        """Execute spectrochemical series query."""
        qt = query_type.lower().strip()

        if qt == "full_series":
            return self._query_full_series(geometry)
        elif qt == "compare_ligands":
            return self._compare_ligands(ligands or [], geometry)
        elif qt == "metal_trend":
            return self._query_metal_trend(metal_ion)
        elif qt == "ligand_position":
            return self._ligand_positions(ligands or [])
        elif qt == "delta_estimate":
            return self._estimate_delta(metal_ion, ligands, geometry)
        else:
            raise ChemMCPError(
                f"Unknown query_type: '{query_type}'. "
                f"Use: 'full_series', 'compare_ligands', 'metal_trend', 'ligand_position', or 'delta_estimate'."
            )

    def _query_full_series(self, geometry: str) -> dict:
        """Return complete spectrochemical series."""
        lig_str = " < ".join([f"{sym}" for _, sym, _, _ in self._ligand_series])
        metal_str = " < ".join([f"{m}" for m, _, _ in self._metal_trend])

        delta_note = ""
        if geometry == "tetrahedral":
            delta_note = "\nNote: For tetrahedral geometry, Δt ≈ (4/9) × Δo. The ordering remains the same."

        explanation = (
            f"The spectrochemical series orders ligands by the magnitude of crystal field splitting they produce.\n"
            f"Weak-field ligands (left) produce small Δ → favor high-spin configurations.\n"
            f"Strong-field ligands (right) produce large Δ → can force low-spin configurations.{delta_note}\n\n"
            f"Key factors affecting Δ:\n"
            f"• Ligand nature (σ-donor/π-acceptor ability)\n"
            f"• Metal oxidation state (higher OS → larger Δ)\n"
            f"• Metal period (5d > 4d > 3d for same ligand)"
        )

        return {
            "query_type": "full_series",
            "series_data": {
                "ligand_series": lig_str,
                "metal_series": metal_str,
                "period_trend": "3d (first row transition metals) < 4d (second row) < 5d (third row)",
                "oxidation_state_rule": "Higher oxidation state → larger Δo (e.g., Δo(Fe³⁺) > Δo(Fe²⁺) by ~30-50%)",
                "delta_t_ratio": "Δt ≈ (4/9) × Δo for same metal + ligands in tetrahedral vs octahedral geometry",
                "ligand_details": [
                    {"name": n, "symbol": s, "relative_factor": f, "category": c}
                    for n, s, f, c in self._ligand_series
                ],
                "metal_details": [
                    {"ion": m, "relative_factor": f, "note": n}
                    for m, f, n in self._metal_trend
                ],
            },
            "explanation": explanation,
            "references": [
                "Figgis, B.N. & Hitchman, M.A. (2000). Ligand Field Theory and Its Applications.",
                "Huheey, J.E. et al. (1993). Inorganic Chemistry, 4th ed.",
                "Shriver & Atkins (2010). Inorganic Chemistry, 5th ed.",
                "Jørgensen, C.K. (1962). Absorption spectra and chemical bonding in complexes.",
            ],
        }

    def _compare_ligands(self, ligands: list, geometry: str) -> dict:
        """Compare specific ligands' field strengths."""
        if not ligands:
            raise ChemMCPError("Please provide at least one ligand to compare.")

        # Find each ligand in the series
        found = []
        not_found = []
        for lg in ligands:
            lg_lower = lg.lower()
            match = None
            for name, symbol, factor, cat in self._ligand_series:
                if lg_lower in name.lower() or lg_lower in symbol.lower() or lg_lower.replace("-", "") == symbol.replace("⁻", "").replace("²", ""):
                    match = {"name": name, "symbol": symbol, "factor": factor, "category": cat}
                    break
            if match:
                found.append(match)
            else:
                not_found.append(lg)

        # Sort found by factor
        found.sort(key=lambda x: x["factor"])

        ordering = " < ".join([f["symbol"] for f in found])
        rel_dict = {f["symbol"]: f["factor"] for f in found}
        cat_dict = {f["symbol"]: f["category"] for f in found}

        explanation = (
            f"Compared {len(found)} ligand(s): {ordering}\n"
            f"Relative to H₂O = 1.00 (reference intermediate field).\n"
            f"Factors < 1.0: weaker than water (weak-field); > 1.0: stronger than water."
        )
        if not_found:
            explanation += f"\nNote: Not found in series: {', '.join(not_found)}"

        if geometry == "tetrahedral":
            explanation += "\nFor tetrahedral: multiply all factors by 4/9 ≈ 0.444."

        return {
            "query_type": "compare_ligands",
            "series_data": {
                "ordering": ordering,
                "relative_delta": rel_dict,
                "field_category": cat_dict,
                "weakest": found[0]["symbol"] if found else None,
                "strongest": found[-1]["symbol"] if found else None,
                "not_found": not_found,
            },
            "explanation": explanation,
            "references": [],
        }

    def _query_metal_trend(self, metal_ion: str = None) -> dict:
        """Return metal ion trend data."""
        if metal_ion:
            m = metal_ion.replace("+", "").upper()
            matches = [entry for entry in self._metal_trend if m in entry[0].replace("+", "").replace("⁺", "")]
            if matches:
                explanation = f"{metal_ion}: relative Δo factor = {matches[0][1]} ({matches[0][2]})"
            else:
                explanation = f"'{metal_ion}' not found in database. General trend: Mn²⁺ < Ni²⁺ < Co²⁺ < Fe²⁺ < V²⁺ < Fe³⁺ < Co³⁺ < Cr³⁺ < 4d/5d metals"
        else:
            explanation = "Metal ions ordered by increasing Δo (with same ligand):"

        return {
            "query_type": "metal_trend",
            "series_data": {
                "ordering": " < ".join([m for m, _, _ in self._metal_trend]),
                "details": [{"ion": m, "factor": f, "note": n} for m, f, n in self._metal_trend],
                "rules": {
                    "oxidation_state": "M³⁺ > M²⁺ > M⁺ (higher charge → larger Δ)",
                    "period": "5d > 4d > 3d (larger radial extension → stronger interaction)",
                    "triad": "Within a group: Δo increases down the group",
                },
            },
            "explanation": explanation,
            "references": [],
        }

    def _ligand_positions(self, ligands: list) -> dict:
        """Find position/rank of specific ligands in the series."""
        positions = {}
        for lg in ligands:
            for i, (name, symbol, factor, cat) in enumerate(self._ligand_series):
                if lg.lower() in name.lower() or lg.lower() in symbol.lower():
                    positions[symbol] = {
                        "rank": i + 1,
                        "total": len(self._ligand_series),
                        "percentile": round((i + 1) / len(self._ligand_series) * 100),
                        "factor": factor,
                        "category": cat,
                    }
                    break

        return {
            "query_type": "ligand_position",
            "series_data": positions,
            "explanation": f"Position of queried ligands in the full spectrochemical series ({len(self._ligand_series)} total ligands).",
            "references": [],
        }

    def _estimate_delta(self, metal_ion: str, ligands: list, geometry: str) -> dict:
        """Rough Δo estimation based on metal + ligand."""
        if not metal_ion or not ligands:
            raise ChemMCPError("Both metal_ion and ligands needed for delta_estimate.")

        # Base Δo from metal
        base_metal = 10000  # cm⁻¹ default
        m_matches = [entry for entry in self._metal_trend if metal_ion.replace("+","").upper() in entry[0].replace("+","").replace("⁺","")]
        if m_matches:
            base_metal = int(base_metal * m_matches[0][1])

        # Adjust for ligand
        for lg in ligands:
            for name, symbol, factor, _ in self._ligand_series:
                if lg.lower() in name.lower() or lg.lower() in symbol.lower():
                    base_metal = int(base_metal * factor)
                    break

        if geometry == "tetrahedral":
            delta_t = int(base_metal * 4 / 9)
            return {
                "query_type": "delta_estimate",
                "series_data": {
                    "estimated_delta_o": base_metal,
                    "estimated_delta_t": delta_t,
                    "metal": metal_ion,
                    "ligands": ligands,
                    "geometry": geometry,
                    "unit": "cm⁻¹ (approximate)",
                    "disclaimer": "Rough estimate only. Use experimental values when available.",
                },
                "explanation": f"Estimated Δo ≈ {base_metal} cm⁻¹ for [{metal_ion}({'+'.join(ligands)})] in octahedral field; Δt ≈ {delta_t} cm⁻¹ in tetrahedral.",
                "references": [],
            }

        return {
            "query_type": "delta_estimate",
            "series_data": {
                "estimated_delta_o": base_metal,
                "metal": metal_ion,
                "ligands": ligands,
                "geometry": geometry,
                "unit": "cm⁻¹ (approximate)",
                "disclaimer": "Rough estimate only.",
            },
            "explanation": f"Estimated Δo ≈ {base_metal} cm⁻¹ for [{metal_ion}({'+'.join(ligands)})].",
            "references": [],
        }

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        if not parts:
            raise ChemMCPError("Format: 'query_type [ligand1 ligand2 ...] [metal_ion]'. Example: 'full_series' or 'compare_ligands NH3 H2O'")
        qt = parts[0]
        ligands = [p for p in parts[1:] if not p[0].isupper() or p in ("NH3", "H2O", "CN-", "Cl-", "OH-", "en", "CO")]
        metal = parts[-1] if parts and parts[-1][0].isupper() and len(parts[-1]) <= 5 and "+" in parts[-1] else None
        if metal and metal in ligands:
            metal = None
        return self._run_base(qt, ligands if ligands else None, metal)
