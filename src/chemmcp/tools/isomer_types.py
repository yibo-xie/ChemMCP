"""
配合物异构体类型分析工具
Analyze possible isomer types for coordination complexes (geometric, optical, ionization, linkage, hydrate).
"""
import logging
import re
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class IsomerTypes(BaseTool):
    """
    分析配合物可能的异构体类型。
    输入化学式或MAmBn型表示，输出可能的几何/光学/电离/键合/水合异构体。
    """
    __version__ = "0.1.0"
    name = "IsomerTypes"
    func_name = "isomer_types"
    description = "Analyze possible isomer types for coordination complexes: geometric (cis/trans, fac/mer), optical, ionization, linkage, and hydrate isomers."
    implementation_description = "Uses combinatorial rules based on coordination geometry and ligand composition. Handles MA2B4, MA3B3, M(AA)2B2, M(AA)3 type formulas and common named complexes."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Coordination Chemistry", "Isomers", "Stereochemistry", "Geometric Isomers", "Optical Isomers"]
    required_envs = []

    code_input_sig = [
        ("complex_formula", "str", "N/A", "Complex formula or notation, e.g., '[Co(NH3)5Cl]SO4', 'MA3B3', '[Pt(NH3)2Cl2]', 'M(AA)2B2'."),
        ("geometry", "str", "octahedral", "Coordination geometry: 'octahedral', 'tetrahedral', or 'square_planar'."),
        ("analyze_all", "bool", "True", "Whether to analyze all isomer types or only geometric."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query string: 'formula [geometry]', e.g., '[Co(NH3)4Cl2] octahedral' or 'MA2B4'."),
    ]

    output_sig = [
        ("complex_formula", "str", "The input formula analyzed."),
        ("geometry", "str", "Assumed geometry."),
        ("isomer_summary", "dict", "Summary of all possible isomer types with counts and descriptions."),
        ("geometric_isomers", "list", "Details of geometric isomers (if any)."),
        ("optical_isomers", "list", "Details of optical isomers/chirality (if any)."),
        ("other_isomers", "list", "Ionization, linkage, hydrate isomers (if applicable)."),
        ("total_isomer_count", "int", "Total number of possible stereoisomers."),
        ("explanation", "str", "Detailed explanation of the analysis."),
    ]

    examples = [
        {
            "code_input": {
                "complex_formula": "[Co(NH3)4Cl2]+",
                "geometry": "octahedral",
                "analyze_all": True,
            },
            "text_input": {
                "query": "[Co(NH3)4Cl2] octahedral"
            },
            "output": {
                "complex_formula": "[Co(NH3)4Cl2]+",
                "geometry": "Octahedral",
                "isomer_summary": {"geometric": 2, "optical": 0, "ionization": 0, "linkage": 0, "hydrate": 0},
                "geometric_isomers": [{"name": "cis-[Co(NH3)4Cl2]+", "description": "Two Cl ligands adjacent (90°)"}, {"name": "trans-[Co(NH3)4Cl2]+", "description": "Two Cl ligands opposite (180°)"}],
                "optical_isomers": [],
                "other_isomers": [],
                "total_isomer_count": 2,
                "explanation": "MA4B2 octahedral → cis/trans geometric isomers. Neither is chiral (has mirror planes).",
            }
        },
        {
            "code_input": {
                "complex_formula": "[Co(en)3]3+",
                "geometry": "octahedral",
                "analyze_all": True,
            },
            "text_input": {
                "query": "[Co(en)3]3+"
            },
            "output": {
                "complex_formula": "[Co(en)3]3+",
                "geometry": "Octahedral",
                "isomer_summary": {"geometric": 1, "optical": 2, "ionization": 0, "linkage": 0, "hydrate": 0},
                "geometric_isomers": [{"name": "Λ-cis-[Co(en)3]3+", "description": "Left-handed propeller"}, {"name": "Δ-cis-[Co(en)3]3+", "description": "Right-handed propeller"}],
                "optical_isomers": [{"name": "Δ enantiomer", "description": "Right-handed helicity, [α]D > 0"}, {"name": "Λ enantiomer", "description": "Left-handed helicity, [α]D < 0"}],
                "other_isomers": [],
                "total_isomer_count": 2,
                "explanation": "Tris-chelate M(AA)3 is optically active (chiral). Δ and Λ enantiomers. No geometric isomers.",
            }
        },
        {
            "code_input": {
                "complex_formula": "[Co(NH3)5SO4]Br",
                "geometry": "octahedral",
                "analyze_all": True,
            },
            "text_input": {
                "query": "[Co(NH3)5SO4]Br"
            },
            "output": {
                "complex_formula": "[Co(NH3)5SO4]Br",
                "geometry": "Octahedral",
                "isomer_summary": {"geometric": 0, "optical": 0, "ionization": 2, "linkage": 0, "hydrate": 0},
                "geometric_isomers": [],
                "optical_isomers": [],
                "other_isomers": [
                    {"type": "ionization isomer", "name": "[Co(NH3)5SO4]Br", "description": "SO4²⁻ inside coordination sphere as ligand, Br⁻ as counterion"},
                    {"type": "ionization isomer", "name": "[Co(NH3)5Br]SO4", "description": "Br⁻ inside coordination sphere as ligand, SO4²⁻ as counterion"},
                ],
                "total_isomer_count": 2,
                "explanation": "Ionization isomers differ in which anion is coordinated vs which is counterion. Give different precipitates with AgNO₃/BaCl₂ tests.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize isomer rule database."""
        # Geometric isomer counts for common formula types
        self._geo_rules = {
            # Octahedral
            "MA6": (1, False),     # Only one arrangement
            "MA5B": (1, False),
            "MA4B2": (2, False),   # cis/trans
            "MA3B3": (2, False),   # fac/mer
            "MA2B2C2": (5, True),  # multiple + some chiral
            "MA2BC4": (3, False),
            "MABC5": (1, False),
            "MA4bc": (1, False),   # bidentate replaces 2 cis positions
            "M(AA)2B2": (2, True),  # cis chiral, trans achiral
            "M(AA)2BC": (3, True),
            "M(AA)3": (1, True),   # optically active (Δ/Λ)
            "M(AA)B4": (1, False),
            "M(A-B)3": (2, True),  # facial and meridional each with 2 enantiomers? No, fac not chiral
            "M(AB)3": (2, True),   # unsym bidentate: fac (chiral pair) + mer (achiral)
            # Square planar
            "sp-MA4": (1, False),
            "sp-MA3B": (1, False),
            "sp-MA2B2": (2, False), # cis/trans (no chirality in sq planar)
            "sp-MABCD": (3, False), # 3 isomers (none chiral)
            # Tetrahedral
            "tet-MA4": (1, False),
            "tet-MA2B2": (1, False), # all adjacent in tetrahedron
            "tet-MABCD": (2, True),  # chiral pair
        }

        # Known complex → formula type mapping
        self._known_complexes = {
            "[co(nh3)4cl2]": ("MA4B2", "octahedral"),
            "[co(nh3)3cl3]": ("MA3B3", "octahedral"),
            "[pt(nh3)2cl2]": ("MA2B2", "square_planar"),
            "[pd(nh3)2cl2]": ("MA2B2", "square_planar"),
            "[co(en)3]": ("M(AA)3", "octahedral"),
            "[cr(en)2cl2]": ("M(AA)2B2", "octahedral"),
            "[cr(nh3)4cl2]": ("MA4B2", "octahedral"),
            "[ni(cn)4]2-": ("MA4", "square_planar"),
            "[ptcl4]2-": ("MA4", "square_planar"),
            "[co(nh3)5so4]br": ("ionization", "octahedral"),
            "[co(nh3)5(no2)cl]": ("linkage", "octahedral"),
            "[cr(h2o)5cl]so4": ("ionization", "octahedral"),
            "[cr(h2o)6]cl3": ("hydrate", "octahedral"),
            "[co(h2o)6]cl3": ("hydrate", "octahedral"),
        }

    def _classify_formula(self, formula: str, geometry: str) -> tuple:
        """Classify a formula into a known pattern type."""
        raw_lower = formula.lower().replace(" ", "")
        # Strip charge BEFORE removing brackets to avoid confusing charge digits with ligand counts
        # Charges appear outside ] like ...]3+ or ...]+
        f_lower = re.sub(r'\]?[+\-]\d*$', '', raw_lower).replace(']', '')  # '...co(en)33]+' -> '...co(en)3'
        f_alt = re.sub(r'\]?\d*[+\-]$', '', raw_lower).replace(']', '')      # fallback

        # Check known complexes first (normalize both sides: strip brackets)
        f_normalized = f_lower.replace("[", "").replace("]", "")
        f_norm_alt = f_alt.replace("[", "").replace("]", "")
        for key, (pattern, geo) in self._known_complexes.items():
            key_norm = key.replace("[", "").replace("]", "")
            if key_norm == f_normalized or key_norm == f_norm_alt:
                return pattern, geo

        # Check notation patterns
        if f_lower.startswith("m("):
            return f_lower, geometry.lower()

        # Try to parse bracketed formulas
        # Count different ligand types
        ligands = re.findall(r'[A-Z][a-z]?[\d]*', f_lower.split("]").pop() if "]" in f_lower else f_lower)

        # Simplified classification by unique ligand count
        unique_ligands = set(re.findall(r'[A-Z][a-z]?|\([a-z]+\)', f_lower))
        n_unique = len(unique_ligands) - 1  # subtract metal

        if geometry.lower() == "square_planar":
            if n_unique <= 1:
                return "sp-MA4", "square_planar"
            elif n_unique == 2:
                return "sp-MA2B2", "square_planar"
            else:
                return "sp-MABCD", "square_planar"

        if geometry.lower() == "tetrahedral":
            if n_unique <= 2:
                return "tet-MA2B2" if n_unique == 2 else "tet-MA4", "tetrahedral"
            elif n_unique >= 4:
                return "tet-MABCD", "tetrahedral"

        # Default octahedral classification
        if n_unique <= 1:
            return "MA6", "octahedral"
        elif n_unique == 2:
            # Could be MA4B2 or MA3B3
            return "MA4B2", "octahedral"  # default assumption
        elif n_unique == 3:
            return "MA2B2C2", "octahedral"
        elif n_unique >= 4:
            return "MABC5", "octahedral"

        return "unknown", geometry.lower()

    def _analyze_geometric(self, pattern: str, geometry: str, formula: str) -> list:
        """Analyze geometric isomers."""
        results = []
        geo = geometry.lower()

        if pattern == "MA4B2" and geo == "octahedral":
            results = [
                {"name": "cis-isomer", "description": "Two B ligands at 90° to each other. Dipole moment ≠ 0."},
                {"name": "trans-isomer", "description": "Two B ligands opposite (180°). Dipole moment may be 0 if symmetric."},
            ]
        elif pattern == "MA3B3" and geo == "octahedral":
            results = [
                {"name": "facial (fac) isomer", "description": "Three B ligands occupying one face of the octahedron (mutually cis, 90°)."},
                {"name": "meridional (mer) isomer", "description": "Three B ligands in a plane through the metal (two trans at 180°, one cis at 90°)."},
            ]
        elif pattern == "M(AA)2B2" and geo == "octahedral":
            results = [
                {"name": "cis-M(AA)2B2", "description": "Two B ligands cis; chiral (exists as Δ/Λ pair)."},
                {"name": "trans-M(AA)2B2", "description": "Two B ligands trans; achiral (mirror plane present)."},
            ]
        elif pattern == "M(AA)3" and geo == "octahedral":
            results = [
                {"name": "Δ-[M(AA)3]", "description": "Right-handed propeller arrangement of chelate rings."},
                {"name": "Λ-[M(AA)3]", "description": "Left-handed propeller arrangement of chelate rings."},
            ]
        elif pattern == "sp-MA2B2":
            results = [
                {"name": "cis-square planar", "description": "Two B ligands adjacent (90°). Usually more reactive."},
                {"name": "trans-square planar", "description": "Two B ligands opposite (180°). More stable typically."},
            ]
        elif pattern == "sp-MABCD":
            results = [
                {"name": "isomer 1", "description": "A trans to B"},
                {"name": "isomer 2", "description": "A trans to C"},
                {"name": "isomer 3", "description": "A trans to D"},
            ]
        elif pattern == "tet-MABCD":
            results = [
                {"name": "Δ enantiomer", "description": "Chiral tetrahedral complex with 4 different ligands."},
                {"name": "Λ enantiomer", "description": "Mirror image of above."},
            ]

        return results

    def _analyze_optical(self, pattern: str, geometry: str) -> list:
        """Check for optical activity (chirality)."""
        chiral_patterns = {
            "M(AA)3": (True, "Tris-chelate has D₃ symmetry, no Sn, no σ → chiral (Δ/Λ)."),
            "M(AA)2B2": (True, "cis isomer lacks improper rotation axis → chiral; trans is achiral."),
            "M(AB)3": (True, "fac isomer with unsymmetrical bidentate is chiral; mer is achiral."),
            "MA2B2C2": (True, "Some arrangements are chiral (all-cis with no pairs trans)."),
            "tet-MABCD": (True, "Tetrahedral with 4 different monodentate ligands → chiral."),
        }

        is_chiral, reason = chiral_patterns.get(pattern, (False, ""))
        if is_chiral:
            return [
                {"name": "Δ (delta) enantiomer", "description": f"Right-handed configuration. {reason}"},
                {"name": "Λ (lambda) enantiomer", "description": f"Left-handed configuration (mirror image of Δ). {reason}"},
            ]
        return []

    def _analyze_other(self, formula: str) -> list:
        """Check for ionization, linkage, hydrate isomers."""
        others = []
        f_lower = formula.lower().replace(" ", "")

        # Ionization isomer indicators: two different anions, one inside/outside
        ionization_patterns = [
            ("[co(nh3)5so4]br", "[Co(NH3)5SO4]Br ↔ [Co(NH3)5Br]SO4"),
            ("[co(nh3)5no2]cl", "[Co(NH3)5(NO2)]Cl ↔ [Co(NH3)5Cl]NO2"),
            ("[pt(nh3)4cl2]br2", "[Pt(NH3)4Cl2]Br2 ↔ [Pt(NH3)4Br2]Cl2 ↔ [Pt(NH3)4ClBr]ClBr"),
            ("[cr(h2o)5cl]so4", "[Cr(H2O)5Cl]SO4 ↔ [Cr(H2O)5SO4]Cl"),
        ]

        for pattern, desc in ionization_patterns:
            if pattern in f_lower:
                names = desc.split(" ↔ ")
                others.append({"type": "ionization isomer", "names": names})
                break

        # Linkage isomer indicators: ambidentate ligands (NO2, SCN, CN, etc.)
        linkage_indicators = ["no2", "scn", "cn", "seo", "oco"]
        for ind in linkage_indicators:
            if ind in f_lower:
                if ind == "no2":
                    others.append({
                        "type": "linkage isomer",
                        "names": ["nitro isomer (M-NO2, N-bonded)", "nitrito isomer (M-ONO, O-bonded)"],
                        "example": "[Co(NH3)5(NO2)]Cl2 ↔ [Co(NH3)5(ONO)]Cl2"
                    })
                elif ind == "scn":
                    others.append({
                        "type": "linkage isomer",
                        "names": ["thiocyanato-S (M-SCN)", "thiocyanato-N (M-NCS)"],
                        "example": "[Co(NH3)5(SCN)]Cl2 ↔ [Co(NH3)5(NCS)]Cl2"
                    })
                break

        # Hydrate isomer indicators
        if "h2o" in f_lower and ("cl" in f_lower or "br" in f_lower or "so4" in f_lower):
            others.append({
                "type": "hydrate (or solvate) isomer",
                "names": ["Water inside coordination sphere as aqua ligand", "Water outside as water of crystallization"],
                "example": "[Cr(H2O)6]Cl3 ↔ [Cr(H2O)5Cl]·H2O·Cl2 ↔ [Cr(H2O)4Cl2]·2H2O·Cl"
            })

        # Coordination isomer (for salts with both cationic and anionic complexes)
        if "[" in formula and formula.count("[") >= 1:
            # Check if it could be a salt of two complexes
            parts = formula.replace("]", "[").split("[")
            if len(parts) > 2:
                others.append({
                    "type": "coordination isomer (possible)",
                    "names": ["Distribution of ligands between cationic and anionic complexes differs"],
                })

        return others

    def _run_base(self, complex_formula: str, geometry: str = "octahedral",
                  analyze_all: bool = True) -> dict:
        """Analyze isomer types for a coordination complex."""
        geo = geometry.lower()
        pattern, detected_geo = self._classify_formula(complex_formula, geo)

        # Geometric isomers
        geo_isomers = self._analyze_geometric(pattern, detected_geo, complex_formula)
        geo_count = len(geo_isomers)

        # Optical isomers
        opt_isomers = self._analyze_optical(pattern, detected_geo) if analyze_all else []
        opt_count = len(opt_isomers)

        # Other isomer types
        other_isomers = self._analyze_other(complex_formula) if analyze_all else []

        # Total count
        total = geo_count + opt_count

        # Build summary
        summary = {
            "geometric": geo_count,
            "optical": opt_count,
            "ionization": sum(1 for o in other_isomers if o["type"] == "ionization isomer"),
            "linkage": sum(1 for o in other_isomers if o["type"] == "linkage isomer"),
            "hydrate": sum(1 for o in other_isomers if "hydrate" in o["type"]),
        }

        # Explanation
        explanation_parts = [
            f"Formula: {complex_formula}",
            f"Detected pattern: {pattern} ({detected_geo.capitalize()} geometry)",
            f"",
        ]
        if geo_count > 0:
            explanation_parts.append(f"Geometric isomers: {geo_count} — {', '.join(i['name'] for i in geo_isomers)}")
        if opt_count > 0:
            explanation_parts.append(f"Optical isomers (enantiomers): {opt_count} — complex is chiral")
        if other_isomers:
            for o in other_isomers:
                explanation_parts.append(f"{o['type'].capitalize()}: possible (see details)")
        if total == 0:
            explanation_parts.append("No isomers expected for this composition.")

        logger.info(f"Isomer analysis: {complex_formula} ({pattern}) → {total} isomers")

        return {
            "complex_formula": complex_formula,
            "geometry": detected_geo.capitalize(),
            "isomer_summary": summary,
            "geometric_isomers": geo_isomers,
            "optical_isomers": opt_isomers,
            "other_isomers": other_isomers,
            "total_isomer_count": total,
            "explanation": "\n".join(explanation_parts),
        }

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        if not parts:
            raise ChemMCPError("Format: 'formula [geometry]'. Example: '[Co(NH3)4Cl2] octahedral' or 'MA3B3'")
        formula = parts[0]
        geo = parts[1] if len(parts) > 1 else "octahedral"
        return self._run_base(formula, geo)
