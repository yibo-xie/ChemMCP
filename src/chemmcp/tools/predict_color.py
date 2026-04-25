"""
基于d-d跃迁预测配合物颜色工具
Predict coordination complex color based on d-d transition (complementary color theory).
"""
import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PredictColor(BaseTool):
    """
    基于晶体场分裂能（Δ）和互补色原理预测配合物颜色。
    输入金属离子、d电子数、场强等，输出预测颜色和吸收/发射波长。
    """
    __version__ = "0.1.0"
    name = "PredictColor"
    func_name = "predict_color"
    description = "Predict the color of a coordination complex based on d-d transitions using complementary color theory. Estimates λmax from crystal field splitting energy."
    implementation_description = "Uses complementary color wheel: the observed color is complementary to the absorbed light wavelength. Δ (cm⁻¹) → λ (nm) = 10⁷/Δ. Includes known colors of common complexes as reference data."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Coordination Chemistry", "Color", "d-d Transition", "Spectroscopy", "CFT"]
    required_envs = []

    code_input_sig = [
        ("metal_ion", "str", "N/A", "Metal ion, e.g., 'Ti3+', 'Cu2+', 'Ni2+', 'Co2+'."),
        ("d_electron_count", "int", "N/A", "Number of d electrons (e.g., 1 for Ti3+, 9 for Cu2+). If omitted, auto-detected from metal_ion."),
        ("geometry", "str", "octahedral", "'octahedral' or 'tetrahedral'."),
        ("field_strength", "str", "intermediate", "'strong', 'weak', or 'intermediate'. Affects Δ estimation."),
        ("ligands", "str", "H2O", "Ligand(s) for Δo estimation (optional)."),
        ("observed_lambda_max_nm", "float", "None", "Experimentally observed λmax in nm. If provided, overrides estimated value."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query string: 'metal_ion [d_count] geometry field_strength [ligands] [lambda_max]', e.g., 'Ti3+ 1 octahedral weak H2O' or 'Cu2+ 9 octahedral weak H2O 800'."),
    ]

    output_sig = [
        ("metal_ion", "str", "Metal ion analyzed."),
        ("d_count", "int", "d-electron count."),
        ("absorption", "dict", "Absorbed light: wavelength (nm), wavenumber (cm⁻¹), color name, region (UV/vis/IR)."),
        ("observed_color", "dict", "Observed (complementary) color: name, hex code, RGB, description."),
        ("transition_type", "str", "Type of d-d transition (e.g., 'T2g→Eg')."),
        ("known_color_reference", "str", "Known color of similar complex for comparison (if available)."),
        ("explanation", "str", "Full explanation of the prediction methodology and reasoning."),
    ]

    examples = [
        {
            "code_input": {
                "metal_ion": "Ti3+",
                "d_electron_count": 1,
                "geometry": "octahedral",
                "field_strength": "weak",
                "ligands": "H2O",
                "observed_lambda_max_nm": None,
            },
            "text_input": {
                "query": "Ti3+ 1 octahedral weak H2O"
            },
            "output": {
                "metal_ion": "Ti3+",
                "d_count": 1,
                "absorption": {"wavelength_nm": 495, "wavenumber_cm": 20200, "color": "green", "region": "visible"},
                "observed_color": {"name": "purple/violet", "hex": "#7F00FF", "description": "Complementary to green absorption"},
                "transition_type": "²T2g → ²Eg (single d electron)",
                "known_color_reference": "[Ti(H2O)6]3+ is violet-purple",
                "explanation": "Ti³⁺ is d¹. Single transition t2g→eg absorbs green (~500 nm) → appears purple/violet.",
            }
        },
        {
            "code_input": {
                "metal_ion": "Cu2+",
                "d_electron_count": 9,
                "geometry": "octahedral",
                "field_strength": "weak",
                "ligands": "H2O",
                "observed_lambda_max_nm": None,
            },
            "text_input": {
                "query": "Cu2+ 9 octahedral weak H2O"
            },
            "output": {
                "metal_ion": "Cu2+",
                "d_count": 9,
                "absorption": {"wavelength_nm": 810, "wavenumber_cm": 12350, "color": "red/NIR", "region": "red/near-IR"},
                "observed_color": {"name": "blue", "hex": "#007FFF", "description": "Complementary to red-orange absorption"},
                "transition_type": "²Eg → ²T2g (hole transition, d⁹ ≡ d⁻¹)",
                "known_color_reference": "[Cu(H2O)6]2+ is blue",
                "explanation": "Cu²⁺ is d⁹ (one hole in eg). Jahn-Teller distorted octahedron. Broad absorption in red-NIR → blue.",
            }
        },
        {
            "code_input": {
                "metal_ion": "Co2+",
                "d_electron_count": 7,
                "geometry": "tetrahedral",
                "field_strength": "weak",
                "ligands": "Cl-",
                "observed_lambda_max_nm": None,
            },
            "text_input": {
                "query": "Co2+ 7 tetrahedral weak Cl-"
            },
            "output": {
                "metal_ion": "Co2+",
                "d_count": 7,
                "absorption": {"wavelength_nm": 680, "wavenumber_cm": 14700, "color": "red", "region": "visible"},
                "observed_color": {"name": "blue-green / teal", "hex": "#008080", "description": "Complementary to red absorption"},
                "transition_type": "⁴A₂ → ⁴T1(P) (tetrahedral d⁷)",
                "known_color_reference": "[CoCl4]2- is deep blue",
                "explanation": "Co²⁺ d⁷ tetrahedral often has intense color due to lack of center of symmetry (relaxing Laporte rule partially).",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize color database."""
        # d-electron counts for common ions
        self._d_counts = {
            "ti3+": 1, "v3+": 2, "cr3+": 3, "cr2+": 4, "mn2+": 5,
            "fe3+": 5, "fe2+": 6, "co3+": 6, "co2+": 7, "ni2+": 8,
            "cu2+": 9, "cu+": 10, "zn2+": 10, "mn3+": 4, "v2+": 3,
        }

        # Known Δo values (cm⁻¹) for common complexes → λmax = 10^7/Δ
        self._delta_data = {
            ("ti3+", "h2o"): 20400,   # [Ti(H2O)6]3+: ~490 nm, purple
            ("v3+", "h2o"): 17800,    # [V(H2O)6]3+: ~560 nm, yellow-green
            ("cr3+", "h2o"): 17400,   # [Cr(H2O)6]3+: ~575 nm, violet
            ("cr3+", "nh3"): 21500,   # [Cr(NH3)6]3+: ~465 nm, yellow
            ("co2+", "h2o"): 9300,    # [Co(H2O)6]2+: ~1075 nm, pink (faint)
            ("co2+", "cl"): 3300,     # [CoCl4]2- tet: ~3000 cm⁻¹? No, use ~600nm
            ("ni2+", "h2o"): 8500,    # [Ni(H2O)6]2+: ~1176 nm, green
            ("cu2+", "h2o"): 12600,   # [Cu(H2O)6]2+: ~794 nm, blue
            ("fe3+", "h2o"): 13700,   # very pale violet
            ("mn2+", "h2o"): 7500,    # very pale pink
        }

        # Complementary color wheel: absorbed wavelength (nm) → observed color
        # Based on artist's color wheel with approximate ranges
        self._complementary = [
            # (absorbed_lambda_min, absorbed_lambda_max, abs_color_name, obs_name, obs_hex)
            (380, 450, "violet", "Yellow-green", "#BFFF00"),
            (450, 485, "blue", "Orange", "#FF7F00"),
            (485, 505, "cyan-green", "Purple/Violet", "#9932CC"),
            (505, 550, "green", "Purple/Magenta", "#FF00FF"),
            (550, 585, "yellow-green", "Violet", "#7F00FF"),
            (585, 620, "yellow", "Blue", "#007FFF"),
            (620, 660, "orange", "Blue-Green", "#00A0FF"),
            (660, 750, "red", "Cyan/Green", "#00FFFF"),
            (750, 1000, "near-IR", "Blue-Green (pale)", "#80C0FF"),
            (1000, 2000, "IR", "Colorless/pale", "#E0E0E0"),
        ]

        # Known colors of common complexes (reference database)
        self._known_colors = {
            "[Ti(H2O)6]3+": ("purple-violet", "#9932CC"),
            "[V(H2O)6]3+": ("green-yellow", "#C8Dc33"),
            "[Cr(H2O)6]3+": ("violet", "#8A2BE2"),
            "[Cr(NH3)6]3+": ("yellow", "#FFD700"),
            "[Cr(CN)6]3-": ("colorless/near-colorless", "#F0F0F0"),
            "[Mn(H2O)6]2+": ("very pale pink", "#FFE4E1"),
            "[Fe(H2O)6]2+": ("very pale green", "#E0FFE0"),
            "[Fe(H2O)6]3+": ("very pale violet/lilac", "#E6E6FA"),
            "[Co(H2O)6]2+": ("pink", "#FFB6C1"),
            "[Co(NH3)6]2+": ("tan/yellow-brown", "#D2B48C"),
            "[CoCl4]2-": ("deep blue", "#00008B"),
            "[Ni(H2O)6]2+": ("green", "#228B22"),
            "[Ni(NH3)6]2+": ("blue-purple", "#6A5ACD"),
            "[Cu(H2O)6]2+": ("blue", "#007FFF"),
            "[Cu(NH3)4]2+": ("deep blue", "#0000CD"),
            "[Zn(H2O)6]2+": ("colorless", "#FFFFFF"),
            "[Fe(CN)6]4-": ("pale yellow", "#FAFAD2"),
            "[Fe(CN)6]3-": ("deep red/orange", "#FF4500"),
            "[Co(NH3)6]3+": ("orange/yellow", "#FFA500"),
            "[PtCl4]2-": ("colorless", "#FFFFFF"),
            "[PdCl4]2-": ("colorless", "#FFFFFF"),
        }

    def _get_d_count(self, metal: str, provided: int = None) -> int:
        if provided is not None:
            return provided
        key = metal.lower().replace(" ", "")
        if key in self._d_counts:
            return self._d_counts[key]
        raise ChemMCPError(f"Unknown metal '{metal}' and no d_count provided.")

    def _estimate_lambda_max(self, metal: str, ligands: str, d_n: int, geo: str) -> int:
        """Estimate λmax in nm from Δo data or generic estimate."""
        key = (metal.lower(), ligands.lower())
        if key in self._delta_data:
            delta = self._delta_data[key]
        else:
            ref_key = (metal.lower(), "h2o")
            if ref_key in self._delta_data:
                delta = self._delta_data[ref_key]
            else:
                # Generic estimate based on d-count and period
                delta = 12000 + d_n * 500  # rough

        if geo == "tetrahedral":
            delta = int(delta * 4 / 9)

        # λ(nm) = 10^7 / Δ(cm⁻¹)
        lam = int(10_000_000 / delta)
        return max(lam, 200)  # sanity: not below UV

    def _complementary_color(self, lambda_nm: int) -> dict:
        """Find complementary color for given absorbed wavelength."""
        for lam_min, lam_max, abs_name, obs_name, obs_hex in self._complementary:
            if lam_min <= lambda_nm < lam_max:
                return {
                    "name": obs_name,
                    "hex_code": obs_hex,
                    "absorbed_range": f"{lam_min}-{lam_max} nm",
                    "absorbed_color": abs_name,
                    "absorbed_lambda_nm": lambda_nm,
                }

        if lambda_nm >= 1000:
            return {"name": "colorless or very pale", "hex_code": "#F5F5F5"}
        return {"name": "unknown (out of visible range)", "hex_code": "#CCCCCC"}

    def _transition_label(self, d_n: int, geo: str) -> str:
        """Generate transition label."""
        if geo == "tetrahedral":
            terms = {1:"⁴T2→⁴A2", 2:"⁴T1(F)→⁴A2", 3:"⁴T1(F)→⁴A2", 4:"⁴A2→⁴T1",
                     5:"⁴T1→⁴T2", 6:"⁴T1→⁴A2", 7:"⁴A2→⁴T1(P)", 8:"⁴T1→⁴A2"}
            return terms.get(d_n, f"d^{d_n} tetrahedral transition")

        oct_terms = {
            1: "²T2g → ²Eg (d¹)",
            2: "³T1g → ³T2g (d²)",
            3: "⁴A2g → ⁴T2g (d³)",
            4: "⁵Eg → ⁵T2g (HS d⁴) / ³T1g → ³T2g (LS d⁴)",
            5: "⁶A1g → ⁶T1g (HS d⁵) / ²T2g → ²Eg (LS d⁵)",
            6: "⁵T2g → ⁵Eg (HS d⁶) / ¹A1g → ¹T1g (LS d⁶)",
            7: "⁴T1g → ⁴T2g (HS d⁷) / ²Eg → ²T2g (LS d⁷)",
            8: "³A2g → ³T2g (d⁸)",
            9: "²Eg → ²T2g (d⁹, Jahn-Teller active)",
            10: "S₀ → S₁ (LMCT/MLCT, d¹⁰ no d-d)",
        }
        return oct_terms.get(d_n, f"d^{d_n} octahedral transition")

    def _run_base(self, metal_ion: str, d_electron_count: int = None,
                  geometry: str = "octahedral", field_strength: str = "intermediate",
                  ligands: str = "H2O", observed_lambda_max_nm: float = None) -> dict:
        """Predict complex color."""
        geo = geometry.lower()
        d_n = self._get_d_count(metal_ion, d_electron_count)

        # Determine λmax
        if observed_lambda_max_nm is not None:
            lam = int(observed_lambda_max_nm)
            delta = int(10_000_000 / lam) if lam > 0 else 0
        else:
            lam = self._estimate_lambda_max(metal_ion, ligands, d_n, geo)
            delta = int(10_000_000 / lam)

        # Complementary color
        comp = self._complementary_color(lam)

        # Spectral region
        if lam < 380:
            region = "UV (ultraviolet)"
        elif lam <= 750:
            region = "visible"
        elif lam <= 2500:
            region = "near-IR (infrared)"
        else:
            region = "IR"

        # Absorption color name
        abs_color = self._absorb_color_name(lam)

        # Look up known complex color
        complex_formula = f"[{metal_ion}({ligands})]"
        known = self._known_colors.get(complex_formula)
        known_str = f"{complex_formula} is known to be {known[0]} ({known[1]})" if known else "No exact reference found in database"

        # Transition type
        trans = self._transition_label(d_n, geo)

        # Explanation
        explanation = (
            f"{metal_ion} is d^{d_n}. "
            f"In {geo} geometry with {ligands}: "
            f"estimated Δ ≈ {delta} cm⁻¹ → λmax ≈ {lam} nm ({region}).\n"
            f"This wavelength falls in the {abs_color} region of the spectrum.\n"
            f"The complementary (observed) color is {comp['name']}.\n\n"
            f"Transition: {trans}\n"
            f"{known_str}"
        )

        logger.info(f"Color prediction: {metal_ion} d^{d_n} → absorbs {lam}nm ({abs_color}), appears {comp['name']}")

        return {
            "metal_ion": metal_ion,
            "d_count": d_n,
            "absorption": {
                "wavelength_nm": lam,
                "wavenumber_cm": delta,
                "color": abs_color,
                "region": region,
            },
            "observed_color": {
                "name": comp["name"],
                "hex_code": comp["hex_code"],
                "description": f"Complementary to absorbed {abs_color} light ({lam} nm)",
            },
            "transition_type": trans,
            "known_color_reference": known_str,
            "explanation": explanation,
        }

    def _absorb_color_name(self, lam: int) -> str:
        """Get color name for an absorbed wavelength."""
        color_ranges = [
            (380, 430, "violet"), (430, 455, "blue-violet"), (455, 490, "blue"),
            (490, 520, "cyan/green"), (520, 565, "green"), (565, 590, "yellow-green"),
            (590, 625, "yellow/orange"), (625, 670, "orange-red"), (670, 750, "red"),
            (750, 1000, "near-infrared"), (1000, 9999, "infrared"),
        ]
        for lo, hi, name in color_ranges:
            if lo <= lam < hi:
                return name
        return "far-IR"

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        if len(parts) < 1:
            raise ChemMCPError("Format: 'metal_ion [d_count] geometry [field_strength] [ligands] [lambda_max]'")
        metal = parts[0]
        d_count = None
        geo = "octahedral"
        field = "intermediate"
        ligands = "H2O"
        lam = None

        i = 1
        while i < len(parts):
            p = parts[i]
            if p.isdigit() and d_count is None and int(p) <= 10:
                d_count = int(p)
            elif p in ("octahedral", "tetrahedral"):
                geo = p
            elif p in ("strong", "weak", "intermediate"):
                field = p
            elif lam is None and self._is_float(p):
                lam = float(p)
            else:
                ligands = p
            i += 1

        kwargs = {"metal_ion": metal, "geometry": geo, "field_strength": field, "ligands": ligands}
        if d_count is not None:
            kwargs["d_electron_count"] = d_count
        if lam is not None:
            kwargs["observed_lambda_max_nm"] = lam
        return self._run_base(**kwargs)

    @staticmethod
    def _is_float(s: str) -> bool:
        try:
            float(s)
            return True
        except ValueError:
            return False
