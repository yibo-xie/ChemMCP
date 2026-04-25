"""
晶体场分裂能分析工具 (Δo / Δt)
Crystal field splitting energy analysis for coordination complexes.
"""
import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CrystalFieldSplitting(BaseTool):
    """
    分析配合物的晶体场分裂能（八面体Δo或四面体Δt）。
    计算CFSE、电子排布、自旋状态等。
    """
    __version__ = "0.1.0"
    name = "CrystalFieldSplitting"
    func_name = "crystal_field_splitting"
    description = "Analyze crystal field splitting energy (Δo for octahedral, Δt for tetrahedral) including CFSE calculation, electron configuration, and spin state."
    implementation_description = "Uses spectrochemical series data and literature Δo values for common complexes. Computes CFSE based on electron occupancy of split d-orbitals. Supports octahedral, tetrahedral, and square planar geometries."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Coordination Chemistry", "Crystal Field Theory", "CFT", "Splitting", "CFSE"]
    required_envs = []

    code_input_sig = [
        ("metal_ion", "str", "N/A", "Metal ion with oxidation state, e.g., 'Cr3+', 'Fe2+', 'Co3+'."),
        ("geometry", "str", "octahedral", "Geometry: 'octahedral', 'tetrahedral', or 'square_planar'."),
        ("ligand_field_strength", "str", "weak", "Ligand field strength: 'strong', 'weak', or 'intermediate'. Affects spin state for d4-d7."),
        ("ligands", "str", "H2O", "Specific ligand(s) for more accurate Δo estimation (optional)."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query string: 'metal_ion geometry field_strength [ligands]', e.g., 'Cr3+ octahedral weak H2O' or 'Co3+ octahedral strong NH3'."),
    ]

    output_sig = [
        ("metal_ion", "str", "Metal ion analyzed."),
        ("d_electron_count", "int", "Number of d electrons."),
        ("geometry", "str", "Coordination geometry."),
        ("splitting_parameter", "dict", "Δ value in cm⁻¹, kJ/mol, and eV with source/explanation."),
        ("electron_configuration", "str", "d-orbital electron configuration (e.g., 't2g³ eg²')."),
        ("cfse", "dict", "Crystal Field Stabilization Energy in Δo units, cm⁻¹, kJ/mol."),
        ("spin_state", "str", "'high-spin' or 'low-spin' (for octahedral d4-d7)."),
        ("unpaired_electrons", "int", "Number of unpaired electrons."),
        ("orbital_diagram", "str", "ASCII/text orbital diagram."),
        ("explanation", "str", "Detailed explanation of the analysis."),
    ]

    examples = [
        {
            "code_input": {
                "metal_ion": "Cr3+",
                "geometry": "octahedral",
                "ligand_field_strength": "weak",
                "ligands": "H2O",
            },
            "text_input": {
                "query": "Cr3+ octahedral weak H2O"
            },
            "output": {
                "metal_ion": "Cr3+",
                "d_electron_count": 3,
                "geometry": "octahedral",
                "splitting_parameter": {"delta_cm": 17400, "delta_kjmol": 208.0, "delta_ev": 2.16, "source": "[Cr(H2O)6]3+ experimental"},
                "electron_configuration": "t2g³ eg⁰",
                "cfse": {"delta_units": "-1.2 Δo", "cm_1": -20880, "kjmol": -249.6},
                "spin_state": "high-spin (only option for d3)",
                "unpaired_electrons": 3,
                "orbital_diagram": "eg:   _   _\nt2g: ↑   ↑   ↑",
                "explanation": "Cr3+ is d3. In octahedral field, all 3 electrons occupy t2g orbitals with parallel spins (Hund's rule). Maximum stabilization.",
            }
        },
        {
            "code_input": {
                "metal_ion": "Fe2+",
                "geometry": "octahedral",
                "ligand_field_strength": "strong",
                "ligands": "CN-",
            },
            "text_input": {
                "query": "Fe2+ octahedral strong CN-"
            },
            "output": {
                "metal_ion": "Fe2+",
                "d_electron_count": 6,
                "geometry": "octahedral",
                "splitting_parameter": {"delta_cm": 33000, "delta_kjmol": 394.5, "delta_ev": 4.09, "source": "[Fe(CN)6]4- estimated"},
                "electron_configuration": "t2g⁶ eg⁰",
                "cfse": {"delta_units": "-2.4 Δo", "cm_1": -79200, "kjmol": -946.8},
                "spin_state": "low-spin",
                "unpaired_electrons": 0,
                "orbital_diagram": "eg:   _   _\nt2g: ↓↑  ↓↑  ↓↑",
                "explanation": "Fe2+ is d6 with strong-field ligands (CN-). All 6 electrons pair in t2g orbitals (low-spin). Large CFSE gain compensates pairing energy.",
            }
        },
        {
            "code_input": {
                "metal_ion": "Co2+",
                "geometry": "tetrahedral",
                "ligand_field_strength": "weak",
                "ligands": "Cl-",
            },
            "text_input": {
                "query": "Co2+ tetrahedral weak Cl-"
            },
            "output": {
                "metal_ion": "Co2+",
                "d_electron_count": 7,
                "geometry": "tetrahedral",
                "splitting_parameter": {"delta_cm": 3300, "delta_kjmol": 39.4, "delta_ev": 0.41, "source": "Δt ≈ 4/9 Δo(oct), estimated"},
                "electron_configuration": "e⁴ t2³",
                "cfse": {"delta_units": "-0.6 Δt", "cm_1": -1980, "kjmol": -23.7},
                "spin_state": "always high-spin (tetrahedral)",
                "unpaired_electrons": 3,
                "orbital_diagram": "t2: ↑   ↑   ↑\ne :  ↓↑  ↓↑",
                "explanation": "Co2+ is d7 in tetrahedral field. Δt is small (~4/9 Δo), so always high-spin. e orbitals fill first (lower energy in tetrahedral).",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize CFT database."""
        # d-electron count for common metal ions
        self._d_counts = {
            # Sc(III) → d0, Ti(IV) → d0
            "ti3+": 1, "v3+": 2, "cr3+": 3, "cr2+": 4,
            "mn3+": 4, "mn2+": 5, "fe3+": 5, "fe2+": 6,
            "co3+": 6, "co2+": 7, "ni2+": 8, "ni3+": 7,
            "cu2+": 9, "cu+": 10, "zn2+": 10, "ag+": 10,
            "pt2+": 8, "pd2+": 8, "au3+": 8,
            "ti2+": 2, "v2+": 3, "cr+": 5,
            "mo3+": 3, "ru3+": 5, "rh3+": 6,
            "os3+": 5, "ir3+": 6,
        }

        # Experimental Δo values (in cm⁻¹) for selected complexes
        # Format: (metal_ion, ligand_key) → delta_o (cm^-1)
        self._delta_o_data = {
            # Cr(III) complexes
            ("cr3+", "h2o"): 17400, ("cr3+", "nh3"): 21500,
            ("cr3+", "en"): 21900, ("cr3+", "cn-"): 26600,
            ("cr3+", "f-"): 14900, ("cr3+", "cl-"): 13800,
            ("cr3+", "oh-"): 17000,
            # Co(III) complexes
            ("co3+", "f-"): 13000, ("co3+", "h2o"): 18600,
            ("co3+", "nh3"): 23000, ("co3+", "en"): 23400,
            ("co3+", "cn-"): 34000,
            # Co(II) complexes
            ("co2+", "h2o"): 9300, ("co2+", "nh3"): 10100,
            # Rh(III)
            ("rh3+", "h2o"): 27000, ("rh3+", "nh3"): 33900,
            # Ir(III)
            ("ir3+", "h2o"): 36000, ("ir3+", "nh3"): 41000,
            # Fe(II/III)
            ("fe2+", "h2o"): 10400, ("fe2+", "nh3"): 10700,
            ("fe2+", "cn-"): 33000, ("fe2+", "h2o_hs"): 10400,
            ("fe3+", "h2o"): 13700, ("fe3+", "cn-"): 35000,
            ("fe3+", "f-"): 15000,
            # Ni(II)
            ("ni2+", "h2o"): 8500, ("ni2+", "nh3"): 10800,
            ("ni2+", "en"): 11600,
            # Mn(II/III)
            ("mn2+", "h2o"): 7500, ("mn3+", "h2o"): 21000,
            # V(II/III)
            ("v2+", "h2o"): 11800, ("v3+", "h2o"): 17700,
            # Cu(II)
            ("cu2+", "h2o"): 12600, ("cu2+", "nh3"): 15100,
            # Pt(II) - square planar reference
            ("pt2+", "cl-"): 26000,  # approximate Δo equivalent
            ("pt2+", "nh3"): 32000,
            # Pd(II)
            ("pd2+", "cl-"): 22400, ("pd2+", "nh3"): 28800,
        }

        # Spectrochemical series positions (relative to H2O = 1.0)
        self._series_factors = {
            "i-": 0.7, "br-": 0.76, "s2-": 0.8, "scn-": 0.76,
            "cl-": 0.80, "no3-": 0.85, "n3-": 0.87, "f-": 0.90,
            "oh-": 0.83, "oxalato": 0.99, "h2o": 1.00,
            "chs-": 1.05, "py": 1.18, "nh3": 1.25,
            "en": 1.28, "bipy": 1.30, "phen": 1.34,
            "no2-": 1.50, "pch3": 1.65, "pr3": 1.70,
            "cn-": 1.70, "co": 1.85,
        }

        # Conversion factors
        self._cm_to_kjmol = 0.01196   # 1 cm⁻¹ = 0.01196 kJ/mol
        self._cm_to_ev = 0.00012398   # 1 cm⁻¹ = 1.23984e-4 eV

    def _get_d_count(self, metal_ion: str) -> int:
        """Get d-electron count from metal ion."""
        key = metal_ion.lower().replace(" ", "")
        if key in self._d_counts:
            return self._d_counts[key]
        raise ChemMCPError(
            f"Unknown metal ion '{metal_ion}'. Known ions: "
            f"{', '.join(sorted(set(k.upper() for k in self._d_counts.keys())))}"
        )

    def _estimate_delta(self, metal_ion: str, ligands: str, geometry: str) -> dict:
        """Estimate Δ value in various units."""
        key = (metal_ion.lower(), ligands.lower())
        factor = self._series_factors.get(ligands.lower(), 1.0)

        if key in self._delta_o_data:
            delta_cm = self._delta_o_data[key]
        else:
            # Try to find same metal with H2O as reference
            ref_key = (metal_ion.lower(), "h2o")
            if ref_key in self._delta_o_data:
                delta_cm = self._delta_o_data[ref_key] * factor
            else:
                # Generic estimate based on metal period and ligand
                delta_cm = 10000 * factor  # rough default ~10000 cm⁻¹
                logger.warning(f"No Δo data for {metal_ion} + {ligands}, using estimate: {delta_cm} cm⁻¹")

        if geometry == "tetrahedral":
            delta_cm = delta_cm * (4.0 / 9.0)  # Δt ≈ 4/9 Δo
        elif geometry == "square_planar":
            # Square planar: Δsp ≈ 1.3–1.5 Δo (very rough)
            delta_cm = int(delta_cm * 1.45)

        return {
            "delta_cm": round(delta_cm),
            "delta_kjmol": round(delta_cm * self._cm_to_kjmol, 1),
            "delta_ev": round(delta_cm * self._cm_to_ev, 2),
            "source": f"{key}" if key in self._delta_o_data else f"estimated from series position ({factor:.2f}× H2O)",
        }

    def _compute_cfse_oct(self, n: int, strong_field: bool) -> dict:
        """
        Compute CFSE for octahedral complex.
        n = number of d electrons.
        Returns cfse in Δo units, orbital config, unpaired count, spin state.
        """
        # Octahedral: t2g (-0.4Δo each), eg (+0.6Δo each)
        # Pairing energy P ≈ 0.8–1.2 Δo (varies by metal); simplified threshold at P ≈ Δo
        configs = {
            # (n, strong?) → (t2g, eg, unpaired, cfse_deltao, spin_label)
            0: (0, 0, 0, 0.0, "N/A"),
            1: (1, 0, 1, -0.4, "N/A"),
            2: (2, 0, 2, -0.8, "N/A"),
            3: (3, 0, 3, -1.2, "N/A"),
            4: (3, 1, 2, -0.6, "high-spin") if not strong_field else (4, 0, 0, -1.6, "low-spin"),
            5: (3, 2, 5, 0.0, "high-spin") if not strong_field else (5, 0, 1, -2.0, "low-spin"),
            6: (4, 2, 4, -0.4, "high-spin") if not strong_field else (6, 0, 0, -2.4, "low-spin"),
            7: (5, 2, 3, -0.8, "high-spin") if not strong_field else (6, 1, 1, -1.8, "low-spin"),
            8: (6, 2, 2, -1.2, "N/A"),  # always same for d8
            9: (6, 3, 1, -0.6, "N/A"),
            10: (6, 4, 0, 0.0, "N/A"),
        }

        if isinstance(configs[4], tuple):
            pass  # already resolved

        entry = configs[n]
        if n in (4, 5, 6, 7):
            entry = configs[n] if isinstance(configs[n], tuple) else (
                configs[n][False] if not strong_field else configs[n][True]
            )

        return {
            "t2g": entry[0], "eg": entry[1],
            "unpaired": entry[2],
            "cfse_deltao": entry[3],
            "spin": entry[4],
        }

    def _compute_cfse_tet(self, n: int) -> dict:
        """
        Compute CFSE for tetrahedral complex (always high-spin).
        Tetrahedral: e (-0.6Δt each), t2 (+0.4Δt each).
        """
        # Always high-spin because Δt is small
        tet_configs = {
            0: (0, 0, 0, 0.0),
            1: (1, 0, 1, -0.6),
            2: (2, 0, 2, -1.2),
            3: (2, 1, 3, -0.8),
            4: (2, 2, 4, -0.4),
            5: (3, 2, 3, -0.6),
            6: (4, 2, 2, -0.8),
            7: (4, 3, 3, -0.6),
            8: (4, 4, 2, -0.4),
            9: (4, 4, 1, -0.2),
            10: (4, 4, 0, 0.0),
        }
        e_count, t2_count, unpaired, cfse = tet_configs[n]
        return {
            "e": e_count, "t2": t2_count,
            "unpaired": unpaired,
            "cfse_deltat": cfse,
        }

    def _make_orbital_diagram_oct(self, t2g: int, eg: int, total_n: int, strong: bool) -> str:
        """Create ASCII orbital diagram for octahedral."""
        # eg orbitals (higher energy)
        eg_str = "   "
        for i in range(min(eg, 2)):
            if i < eg:
                if strong and total_n in (4, 5, 6, 7) and t2g >= 3 + (total_n - 3 if total_n <= 5 else 3):
                    eg_str += "↓↑ "
                else:
                    eg_str += "↑  "
            else:
                eg_str += "_  "

        # t2g orbitals (lower energy)
        t2g_str = ""
        for i in range(3):
            if i < t2g % 3 or (t2g // 3 > i // 3 and t2g % 3 > i % 3):
                if strong and total_n in (5, 6) and t2g > 3:
                    t2g_str += "↓↑ "
                elif t2g > i:
                    t2g_str += "↑  "
                else:
                    t2g_str += "_  "
            elif t2g > 3:
                t2g_str += "↓↑ "
            elif i < t2g:
                t2g_str += "↑  "
            else:
                t2g_str += "_  "

        # Simplified approach
        lines = []
        lines.append("  ┌─────────┐  Energy ↑")
        lines.append("  │  eg     │  ____")
        eg_occ = min(eg, 2)
        occ_str = ""
        for i in range(2):
            if i < eg_occ:
                if strong and total_n in (4, 5, 6, 7) and t2g >= 3:
                    occ_str += "↓↑ "
                else:
                    occ_str += "↑  "
            else:
                occ_str += "_  "
        lines.append(f"  │         │  {occ_str}")
        lines.append("  ├─────────┤  Δo")
        lines.append("  │  t2g    │  ____")
        t2g_occ = min(t2g, 3)
        occ_str2 = ""
        for i in range(3):
            if i < t2g_occ:
                if strong and total_n in (5, 6) and t2g > 3:
                    occ_str2 += "↓↑ "
                else:
                    occ_str2 += "↑  "
            else:
                occ_str2 += "_  "
        lines.append(f"  │         │  {occ_str2}")
        lines.append("  └─────────┘")

        return "\n".join(lines)

    def _run_base(self, metal_ion: str, geometry: str = "octahedral",
                  ligand_field_strength: str = "weak", ligands: str = "H2O") -> dict:
        """Perform complete crystal field splitting analysis."""
        geo = geometry.lower()
        if geo not in ("octahedral", "tetrahedral", "square_planar"):
            raise ChemMCPError(f"Unsupported geometry: '{geometry}'. Use 'octahedral', 'tetrahedral', or 'square_planar'.")

        field = ligand_field_strength.lower()
        if field not in ("strong", "weak", "intermediate"):
            raise ChemMCPError(f"Invalid field strength: '{ligand_field_strength}'. Use 'strong', 'weak', or 'intermediate'.")

        n = self._get_d_count(metal_ion)
        delta_info = self._estimate_delta(metal_ion, ligands, geo)

        if geo == "octahedral":
            is_strong = field == "strong"
            result = self._compute_cfse_oct(n, is_strong)
            t2g, eg = result["t2g"], result["eg"]
            unpaired = result["unpaired"]
            cfse_do = result["cfse_deltao"]
            spin = result["spin"]

            config_str = f"t2g^{t2g} eg^{eg}"
            diagram = self._make_orbital_diagram_oct(t2g, eg, n, is_strong)

            cfse_cm = abs(cfse_do) * delta_info["delta_cm"] if cfse_do < 0 else cfse_do * delta_info["delta_cm"]
            cfse_kj = cfse_cm * self._cm_to_kjmol

            explanation = self._explain_oct(metal_ion, n, is_strong, t2g, eg, unpaired, spin, ligands)

            return {
                "metal_ion": metal_ion,
                "d_electron_count": n,
                "geometry": geo.capitalize(),
                "splitting_parameter": delta_info,
                "electron_configuration": config_str,
                "cfse": {
                    "delta_units": f"{cfse_do:.1f} Δo",
                    "cm_1": round(cfse_do * delta_info["delta_cm"]),
                    "kjmol": round(cfse_do * delta_info["delta_cm"] * self._cm_to_kjmol, 1),
                },
                "spin_state": spin,
                "unpaired_electrons": unpaired,
                "orbital_diagram": diagram,
                "explanation": explanation,
            }

        elif geo == "tetrahedral":
            result = self._compute_cfse_tet(n)
            e_c, t2_c = result["e"], result["t2"]
            unpaired = result["unpaired"]
            cfse_dt = result["cfse_deltat"]

            config_str = f"e^{e_c} t2^{t2_c}"
            diagram = self._make_orbital_diagram_tet(e_c, t2_c)

            explanation = self._explain_tet(metal_ion, n, e_c, t2_c, unpaired, ligands)

            return {
                "metal_ion": metal_ion,
                "d_electron_count": n,
                "geometry": geo.capitalize(),
                "splitting_parameter": delta_info,
                "electron_configuration": config_str,
                "cfse": {
                    "delta_units": f"{cfse_dt:.1f} Δt",
                    "cm_1": round(cfse_dt * delta_info["delta_cm"]),
                    "kjmol": round(cfse_dt * delta_info["delta_cm"] * self._cm_to_kjmol, 1),
                },
                "spin_state": "high-spin (always, for tetrahedral)",
                "unpaired_electrons": unpaired,
                "orbital_diagram": diagram,
                "explanation": explanation,
            }

        else:  # square_planar
            # Simplified square planar (derived from octahedral dz2 removal)
            sp_result = self._compute_cfse_sp(n)
            explanation = self._explain_sp(metal_ion, n, ligands)
            return {
                "metal_ion": metal_ion,
                "d_electron_count": n,
                "geometry": "Square Planar",
                "splitting_parameter": delta_info,
                "electron_configuration": sp_result["config"],
                "cfse": {
                    "delta_units": f"{sp_result['cfse']:.1f} Δ",
                    "cm_1": round(sp_result["cfse"] * delta_info["delta_cm"]),
                    "kjmol": round(sp_result["cfse"] * delta_info["delta_cm"] * self._cm_to_kjmol, 1),
                },
                "spin_state": sp_result.get("spin", "typically diamagnetic"),
                "unpaired_electrons": sp_result.get("unpaired", 0),
                "orbital_diagram": sp_result.get("diagram", "(Square planar diagram - see explanation)"),
                "explanation": explanation,
            }

    def _make_orbital_diagram_tet(self, e_count: int, t2_count: int) -> str:
        """ASCII diagram for tetrahedral (inverted vs octahedral)."""
        lines = []
        lines.append("  ┌─────────┐  Energy ↑")
        lines.append("  │  t2     │  ____")
        occ = ""
        for i in range(3):
            if i < t2_count:
                occ += "↑  "
            else:
                occ += "_  "
        lines.append(f"  │         │  {occ}")
        lines.append("  ├─────────┤  Δt")
        lines.append("  │  e      │  ____")
        occ2 = ""
        for i in range(2):
            if i < e_count:
                occ2 += "↓↑ " if e_count > 2 and i < 2 else "↑  "
            else:
                occ2 += "_  "
        lines.append(f"  │         │  {occ2}")
        lines.append("  └─────────┘")
        return "\n".join(lines)

    def _compute_cfse_sp(self, n: int) -> dict:
        """Simplified square planar CFSE (extreme of octahedral distortion)."""
        # Square planar: dx2-y2 highest, then dxy, then dz2, then dxz/dyz lowest
        # Energy levels (in Δo): dx2-y2: +1.228, dxy: -0.228, dz2: -0.428, dxz=dyz: -0.514
        sp_configs = {
            0: {"config": "d⁰", "cfse": 0.0, "unpaired": 0},
            1: {"config": "(dxz/dyz)¹", "cfse": -1.028, "unpaired": 1},
            2: {"config": "(dxz/dyz)²", "cfse": -2.056, "unpaired": 0},  # usually paired
            3: {"config": "(dxz/dyz)²(dz2)¹", "cfse": -2.484, "unpaired": 1},
            4: {"config": "(dxz/dyz)²(dz2)²", "cfse": -2.912, "unpaired": 0},  # typical Ni(II), Pd(II), Pt(II)
            5: {"config": "(dxz/dyz)²(dz2)²(dxy)¹", "cfse": -3.140, "unpaired": 1},
            6: {"config": "(dxz/dyz)²(dz2)²(dxy)²", "cfse": -3.368, "unpaired": 0},
            7: {"config": "(dxz/dyz)²(dzz)²(dxy)²(dx2-y2)¹", "cfse": -2.140, "unpaired": 1},
            8: {"config": "(dxz/dyz)²(dz2)²(dxy)²(dx2-y2)²", "cfse": -0.912, "unpaired": 0},
            9: {"config": "...(dx2-y2)³", "cfse": +0.316, "unpaired": 1},
            10: {"config": "full", "cfse": 0.0, "unpaired": 0},
        }
        return sp_configs.get(n, {"config": f"d^{n}", "cfse": 0.0, "unpaired": 0})

    def _explain_oct(self, metal, n, strong, t2g, eg, unpaired, spin, ligands):
        """Generate human-readable explanation for octahedral case."""
        field_desc = "strong-field" if strong else "weak-field"
        parts = [
            f"{metal} has {n} d electrons.",
            f"In an octahedral field with {field_desc} ligands ({ligands}), ",
        ]
        if n <= 3 or n >= 8:
            parts.append(f"there is only one way to arrange electrons: t2g^{t2g} eg^{eg}.")
        else:
            hs_or_ls = "low-spin" if strong else "high-spin"
            parts.append(f"the {hs_or_ls} configuration is favored: t2g^{t2g} eg^{eg}.")
        parts.append(f"This gives {unpaired} unpaired electron(s) ({spin}).")
        return "".join(parts)

    def _explain_tet(self, metal, n, e_c, t2_c, unpaired, ligands):
        return (f"{metal} has {n} d electrons. In a tetrahedral field with {ligands}, "
                f"Δt ≈ 4/9 Δo is always small → always high-spin. "
                f"Configuration: e^{e_c} t2^{t2_c}, {unpaired} unpaired electron(s).")

    def _explain_sp(self, metal, n, ligands):
        return (f"{metal} (d{n}) in square planar geometry with {ligands}. "
                f"Square planar is typically favored by d8 metals (Ni2+, Pd2+, Pt2+) with strong-field ligands, "
                f"where large CFSE from dsp2 hybridization overcomes the preference for tetrahedral/octahedral.")

    def _run_text(self, query: str) -> dict:
        """Parse text query."""
        parts = query.strip().split()
        if len(parts) < 2:
            raise ChemMCPError("Format: 'metal_ion geometry [field_strength] [ligands]'. Example: 'Cr3+ octahedral weak H2O'")

        metal_ion = parts[0]
        geometry = parts[1] if len(parts) > 1 else "octahedral"
        field = parts[2] if len(parts) > 2 else "weak"
        ligands = parts[3] if len(parts) > 3 else "H2O"

        return self._run_base(metal_ion, geometry, field, ligands)
