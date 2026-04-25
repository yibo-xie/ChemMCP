"""
配合物磁性预测工具（高自旋/低自旋，磁矩）
Predict magnetism of coordination complexes (high-spin vs low-spin, magnetic moment).
"""
import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PredictMagnetism(BaseTool):
    """
    预测配合物的磁性：高自旋/低自旋、未成对电子数、有效磁矩 μeff。
    使用自旋-only公式 μ = √[n(n+2)] B.M.
    """
    __version__ = "0.1.0"
    name = "PredictMagnetism"
    func_name = "predict_magnetism"
    description = "Predict magnetic properties of coordination complexes: spin state (high/low), number of unpaired electrons, and effective magnetic moment (μeff)."
    implementation_description = "Uses crystal field theory to determine electron configuration, then applies spin-only formula μ = √[n(n+2)] for magnetic moment in Bohr magnetons. Supports octahedral and tetrahedral geometries."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Coordination Chemistry", "Magnetism", "Spin State", "Magnetic Moment", "CFT"]
    required_envs = []

    code_input_sig = [
        ("metal_ion", "str", "N/A", "Metal ion with oxidation state, e.g., 'Fe2+', 'Co3+', 'Mn2+'."),
        ("geometry", "str", "octahedral", "Geometry: 'octahedral' or 'tetrahedral'."),
        ("field_strength", "str", "weak", "Ligand field strength: 'strong' or 'weak' (only matters for octahedral d4-d7)."),
        ("include_orbital_contribution", "bool", "False", "Whether to include orbital angular momentum contribution (rough correction for spin-orbit coupling)."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query string: 'metal_ion geometry field_strength', e.g., 'Fe2+ octahedral strong' or 'Mn2+ tetrahedral weak'."),
    ]

    output_sig = [
        ("metal_ion", "str", "Metal ion analyzed."),
        ("d_electron_count", "int", "Number of d electrons."),
        ("geometry", "str", "Coordination geometry."),
        ("spin_state", "str", "'high-spin', 'low-spin', or 'N/A' (when only one configuration exists)."),
        ("unpaired_electrons", "int", "Number of unpaired electrons."),
        ("spin_only_moment", "float", "Spin-only magnetic moment μso = √[n(n+2)] in Bohr magnetons (B.M.)."),
        ("effective_moment", "float", "Effective magnetic moment μeff (with optional orbital contribution)."),
        ("magnetic_behavior", "str", "'paramagnetic', 'diamagnetic', or 'ferromagnetic' description."),
        ("electron_arrangement", "str", "Description of d-orbital electron arrangement."),
        ("explanation", "str", "Detailed explanation of the prediction."),
    ]

    examples = [
        {
            "code_input": {
                "metal_ion": "Fe2+",
                "geometry": "octahedral",
                "field_strength": "strong",
                "include_orbital_contribution": False,
            },
            "text_input": {
                "query": "Fe2+ octahedral strong"
            },
            "output": {
                "metal_ion": "Fe2+",
                "d_electron_count": 6,
                "geometry": "Octahedral",
                "spin_state": "low-spin",
                "unpaired_electrons": 0,
                "spin_only_moment": 0.0,
                "effective_moment": 0.0,
                "magnetic_behavior": "diamagnetic (low-spin Fe2+, t2g⁶)",
                "electron_arrangement": "(t2g)⁶(eg)⁰ — all electrons paired",
                "explanation": "Fe2+ is d6. With strong-field ligands, Δo > P → all 6 electrons pair in t2g orbitals (t2g⁶). Zero unpaired electrons → diamagnetic.",
            }
        },
        {
            "code_input": {
                "metal_ion": "Fe2+",
                "geometry": "octahedral",
                "field_strength": "weak",
                "include_orbital_contribution": False,
            },
            "text_input": {
                "query": "Fe2+ octahedral weak"
            },
            "output": {
                "metal_ion": "Fe2+",
                "d_electron_count": 6,
                "geometry": "Octahedral",
                "spin_state": "high-spin",
                "unpaired_electrons": 4,
                "spin_only_moment": 4.90,
                "effective_moment": 4.9,
                "magnetic_behavior": "paramagnetic (high-spin Fe2+, t2g⁴ eg²)",
                "electron_arrangement": "(t2g)⁴(eg)² — 4 unpaired electrons",
                "explanation": "Fe2+ is d6. With weak-field ligands, Δo < P → high-spin: t2g⁴ eg² with 4 unpaired electrons. μso = √[4×6] = 4.90 B.M.",
            }
        },
        {
            "code_input": {
                "metal_ion": "Co2+",
                "geometry": "tetrahedral",
                "field_strength": "weak",
                "include_orbital_contribution": False,
            },
            "text_input": {
                "query": "Co2+ tetrahedral"
            },
            "output": {
                "metal_ion": "Co2+",
                "d_electron_count": 7,
                "geometry": "Tetrahedral",
                "spin_state": "high-spin (always for tetrahedral)",
                "unpaired_electrons": 3,
                "spin_only_moment": 3.87,
                "effective_moment": 3.87,
                "magnetic_behavior": "paramagnetic (tetrahedral Co2+, e⁴ t2³)",
                "electron_arrangement": "(e)⁴(t2)³ — 3 unpaired electrons",
                "explanation": "Co2+ is d7. Tetrahedral Δt ≈ 4/9 Δo is always small → always high-spin. e⁴ t2³ with 3 unpaired. μso = √[3×5] = 3.87 B.M.",
            }
        },
        {
            "code_input": {
                "metal_ion": "Mn2+",
                "geometry": "octahedral",
                "field_strength": "weak",
                "include_orbital_contribution": False,
            },
            "text_input": {
                "query": "Mn2+ octahedral weak"
            },
            "output": {
                "metal_ion": "Mn2+",
                "d_electron_count": 5,
                "geometry": "Octahedral",
                "spin_state": "high-spin (only option for d5)",
                "unpaired_electrons": 5,
                "spin_only_moment": 5.92,
                "effective_moment": 5.92,
                "magnetic_behavior": "strongly paramagnetic (high-spin d5, t2g³ eg²)",
                "electron_arrangement": "(t2g)³(eg)² — 5 unpaired electrons (maximum)",
                "explanation": "Mn2+ is d5 (half-filled). Both high-spin and low-spin give same CFSE=0, but pairing energy favors high-spin: t2g³ eg² with 5 unpaired. Maximum paramagnetism. μso = √[5×7] = 5.92 B.M.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize d-electron count database."""
        self._d_counts = {
            "ti3+": 1, "v3+": 2, "cr3+": 3, "cr2+": 4,
            "mn3+": 4, "mn2+": 5, "fe3+": 5, "fe2+": 6,
            "co3+": 6, "co2+": 7, "ni2+": 8, "ni3+": 7,
            "cu2+": 9, "cu+": 10, "zn2+": 10, "ag+": 10,
            "pt2+": 8, "pd2+": 8, "au3+": 8,
            "ti2+": 2, "v2+": 3, "cr+": 5,
            "mo3+": 3, "ru3+": 5, "rh3+": 6,
        }

        # Orbital angular momentum reduction factors (approximate)
        # These account for quenching in different geometries
        self._orbital_factors = {
            "octahedral": {"light": 1.0, "heavy": 1.3},   # light 3d vs heavy 4d/5d
            "tetrahedral": {"light": 1.0, "heavy": 1.2},
        }

    def _get_d_count(self, metal_ion: str) -> int:
        key = metal_ion.lower().replace(" ", "")
        if key not in self._d_counts:
            raise ChemMCPError(
                f"Unknown metal ion '{metal_ion}'. Known: "
                f"{', '.join(sorted(set(k.upper() for k in self._d_counts.keys())))}"
            )
        return self._d_counts[key]

    def _is_heavy_metal(self, metal_ion: str) -> bool:
        """Check if metal is 4d/5d (significant orbital contribution)."""
        heavy = {"pt2+", "pd2+", "au3+", "rh3+", "ir3+", "ru3+", "os3+"}
        return metal_ion.lower() in heavy

    def _determine_spin_oct(self, n: int, strong_field: bool) -> tuple:
        """
        Determine unpaired electrons for octahedral geometry.
        Returns (unpaired, spin_state, t2g, eg).
        """
        # Octahedral electron configurations
        oct_data = {
            0: (0, "N/A", 0, 0),
            1: (1, "N/A", 1, 0),
            2: (2, "N/A", 2, 0),
            3: (3, "N/A", 3, 0),
            4: (2, "low-spin", 4, 0) if strong_field else (4, "high-spin", 3, 1),
            5: (1, "low-spin", 5, 0) if strong_field else (5, "high-spin", 3, 2),
            6: (0, "low-spin", 6, 0) if strong_field else (4, "high-spin", 4, 2),
            7: (1, "low-spin", 6, 1) if strong_field else (3, "high-spin", 5, 2),
            8: (2, "N/A", 6, 2),
            9: (1, "N/A", 6, 3),
            10: (0, "N/A", 6, 4),
        }
        return oct_data[n]

    def _determine_spin_tet(self, n: int) -> tuple:
        """
        Determine unpaired for tetrahedral (always high-spin).
        Returns (unpaired, e_count, t2_count).
        """
        tet_data = {
            0: (0, 0, 0), 1: (1, 1, 0), 2: (2, 2, 0),
            3: (3, 2, 1), 4: (4, 2, 2), 5: (3, 3, 2),
            6: (2, 4, 2), 7: (3, 4, 3), 8: (2, 4, 4),
            9: (1, 4, 4), 10: (0, 4, 4),
        }
        return tet_data[n]

    def _spin_only_moment(self, n_unpaired: int) -> float:
        """μso = √[n(n+2)] in Bohr magnetons."""
        if n_unpaired == 0:
            return 0.0
        return round(math.sqrt(n_unpaired * (n_unpaired + 2)), 2)

    def _run_base(self, metal_ion: str, geometry: str = "octahedral",
                  field_strength: str = "weak",
                  include_orbital_contribution: bool = False) -> dict:
        """Predict magnetic properties."""
        geo = geometry.lower()
        if geo not in ("octahedral", "tetrahedral"):
            raise ChemMCPError("Geometry must be 'octahedral' or 'tetrahedral'.")

        field = field_strength.lower()
        n = self._get_d_count(metal_ion)
        is_heavy = self._is_heavy_metal(metal_ion)

        if geo == "octahedral":
            is_strong = (field == "strong")
            unpaired, spin_state, t2g, eg = self._determine_spin_oct(n, is_strong)
            arr_str = f"(t2g)^{t2g}(eg)^{eg}" if t2g + eg > 0 else "empty"

        else:
            unpaired, e_c, t2_c = self._determine_spin_tet(n)
            spin_state = "high-spin (always, tetrahedral)"
            arr_str = f"(e)^{e_c}(t2)^{t2_c}" if e_c + t2_c > 0 else "empty"

        mu_so = self._spin_only_moment(unpaired)

        # Optional orbital contribution (rough approximation)
        if include_orbital_contribution and unpaired > 0 and is_heavy:
            factor = self._orbital_factors.get(geo, {}).get("heavy", 1.2)
            mu_eff = round(mu_so * factor, 2)
        elif include_orbital_contribution and unpaired > 0:
            factor = self._orbital_factors.get(geo, {}).get("light", 1.0)
            mu_eff = round(mu_so * factor, 2)
        else:
            mu_eff = mu_so

        # Magnetic behavior
        if unpaired == 0:
            behavior = "diamagnetic"
        elif unpaired >= 4:
            behavior = "strongly paramagnetic"
        elif unpaired >= 2:
            behavior = "paramagnetic"
        else:
            behavior = "weakly paramagnetic"

        # Build explanation
        explanation = self._build_explanation(
            metal_ion, n, geo, field, spin_state, unpaired, t2g if geo == "oct" else None,
            eg if geo == "oct" else None, mu_so, mu_eff
        )

        logger.info(f"Magnetism: {metal_ion} ({geo}, {field}) → {spin_state}, {unpaired} unpaired, μ={mu_so} BM")

        return {
            "metal_ion": metal_ion,
            "d_electron_count": n,
            "geometry": geo.capitalize(),
            "spin_state": spin_state,
            "unpaired_electrons": unpaired,
            "spin_only_moment": mu_so,
            "effective_moment": mu_eff,
            "magnetic_behavior": f"{behavior} ({arr_str})",
            "electron_arrangement": f"{arr_str} — {unpaired} unpaired electron(s)",
            "explanation": explanation,
        }

    def _build_explanation(self, metal, n, geo, field, spin, unpaired, t2g=None, eg=None, mu_so=0, mu_eff=0):
        parts = [
            f"{metal} has {n} d electrons.",
            f"In {geo} geometry with {field}-field ligands: ",
        ]
        if n <= 3 or n >= 8:
            parts.append(f"only one electron arrangement exists.")
        else:
            parts.append(f"{spin} configuration is adopted.")
        parts.append(f" {unpaired} unpaired electron(s) → μso = √[{unpaired}×{unpaired+2}] = {mu_so} B.M.")
        if mu_eff != mu_so:
            parts.append(f" With orbital contribution: μeff ≈ {mu_eff} B.M.")
        return "".join(parts)

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        if len(parts) < 1:
            raise ChemMCPError("Format: 'metal_ion [geometry] [field_strength]'. Example: 'Fe2+ octahedral strong'")
        metal = parts[0]
        geo = parts[1] if len(parts) > 1 else "octahedral"
        field = parts[2] if len(parts) > 2 else "weak"
        return self._run_base(metal, geo, field)
