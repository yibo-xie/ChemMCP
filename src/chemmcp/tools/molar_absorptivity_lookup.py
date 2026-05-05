"""
Molar Absorptivity Lookup — 摩尔吸光系数数据库查询
常见化合物在不同波长下的摩尔吸光系数 (ε, M⁻¹cm⁻¹) 数据库
"""
import logging
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ── 摩尔吸光系数数据库 ────────────────────────────────────────────
# 每个条目: {name, mw, data: [{lambda_nm, epsilon, solvent, note}]}
EPSILON_DB: List[dict] = [
    # ═══ Amino Acids & Proteins ═══
    {"name": "Tryptophan", "mw": 204.23, "category": "Amino Acid",
     "data": [
         {"lambda_nm": 280, "epsilon": 5690, "solvent": "H₂O", "note": "π→π* of indole ring; primary protein 280nm contributor."},
         {"lambda_nm": 288, "epsilon": 4750, "solvent": "H₂O", "note": "Secondary peak."},
     ]},
    {"name": "Tyrosine", "mw": 181.19, "category": "Amino Acid",
     "data": [
         {"lambda_nm": 274, "epsilon": 1340, "solvent": "H₂O", "note": "Phenolic π→π*."},
         {"lambda_nm": 222, "epsilon": 8200, "solvent": "H₂O", "note": "Stronger π→π* band."},
         {"lambda_nm": 293, "epsilon": 2300, "solvent": "pH>12 (phenolate)", "note": "Ionized form: greatly enhanced."},
     ]},
    {"name": "Phenylalanine", "mw": 165.19, "category": "Amino Acid",
     "data": [
         {"lambda_nm": 257, "epsilon": 200, "solvent": "H₂O", "note": "Weak benzene B-band."},
         {"lambda_nm": 206, "epsilon": 9300, "solvent": "H₂O", "note": "Strong π→π*."},
     ]},
    {"name": "Cysteine (Cystine)", "mw": 240.48, "category": "Amino Acid",
     "data": [
         {"lambda_nm": 250, "epsilon": 350, "solvent": "H₂O pH 7", "note": "Disulfide bond n→σ*."},
     ]},
    {"name": "Histidine", "mw": 155.16, "category": "Amino Acid",
     "data": [
         {"lambda_nm": 211, "epsilon": 5500, "solvent": "H₂O", "note": "Imidazole π→π*."},
     ]},

    # ═══ Nucleotides & Nucleic Acids ═══
    {"name": "Adenosine / Adenine", "mw": 267.24, "category": "Nucleotide",
     "data": [
         {"lambda_nm": 260, "epsilon": 16400, "solvent": "pH 7 buffer", "note": "Purine base maximum."},
     ]},
    {"name": "Guanosine / Guanine", "mw": 283.24, "category": "Nucleotide",
     "data": [
         {"lambda_nm": 252, "epsilon": 13700, "solvent": "pH 7 buffer", "note": "Purine base."},
         {"lambda_nm": 275, "epsilon": 8000, "solvent": "pH 7 buffer", "note": "Shoulder."},
     ]},
    {"name": "Cytidine / Cytosine", "mw": 243.22, "category": "Nucleotide",
     "data": [
         {"lambda_nm": 267, "epsilon": 6100, "solvent": "pH 7 buffer", "note": "Pyrimidine base."},
     ]},
    {"name": "Thymidine / Thymine", "mw": 242.23, "category": "Nucleotide",
     "data": [
         {"lambda_nm": 264, "epsilon": 7900, "solvent": "pH 7 buffer", "note": "Pyrimidine with methyl."},
     ]},
    {"name": "Uridine / Uracil", "mw": 244.20, "category": "Nucleotide",
     "data": [
         {"lambda_nm": 260, "epsilon": 10200, "solvent": "pH 7 buffer", "note": "Pyrimidine base."},
     ]},
    {"name": "NADH", "mw": 663.43, "category": "Coenzyme",
     "data": [
         {"lambda_nm": 340, "epsilon": 6220, "solvent": "pH 8 buffer", "note": "Reduced nicotinamide ring; key for enzyme assays."},
         {"lambda_nm": 259, "epsilon": 15000, "solvent": "pH 8 buffer", "note": "Adenine moiety."},
     ]},
    {"name": "NAD⁺", "mw": 663.42, "category": "Coenzyme",
     "data": [
         {"lambda_nm": 259, "epsilon": 16900, "solvent": "pH 7 buffer", "note": "Oxidized form; no 340nm absorption."},
     ]},

    # ═══ Common Organic Molecules ═══
    {"name": "Benzene", "mw": 78.11, "category": "Aromatic Hydrocarbon",
     "data": [
         {"lambda_nm": 255, "epsilon": 215, "solvent": "hexane", "note": "B-band (forbidden transition)."},
         {"lambda_nm": 201, "epsilon": 7400, "solvent": "hexane", "note": "E-band (allowed)."},
         {"lambda_nm": 178, "epsilon": 55000, "solvent": "vapor", "note": "Vacuum UV region."},
     ]},
    {"name": "Toluene", "mw": 92.14, "category": "Aromatic Hydrocarbon",
     "data": [
         {"lambda_nm": 261, "epsilon": 300, "solvent": "hexane", "note": "+6nm vs benzene from CH₃ auxochrome."},
         {"lambda_nm": 208, "epsilon": 7900, "solvent": "hexane", "note": "E-band."},
     ]},
    {"name": "Phenol", "mw": 94.11, "category": "Phenol",
     "data": [
         {"lambda_nm": 270, "epsilon": 1450, "solvent": "H₂O", "note": "OH auxochrome enhances B-band."},
         {"lambda_nm": 269, "epsilon": 1450, "solvent": "hexane", "note": "Similar in nonpolar solvent."},
         {"lambda_nm": 287, "epsilon": 2600, "solvent": "NaOH aq", "note": "Phenolate anion: red-shifted + intensified."},
     ]},
    {"name": "Aniline", "mw": 93.13, "category": "Aromatic Amine",
     "data": [
         {"lambda_nm": 230, "epsilon": 8600, "solvent": "H₂O pH 11", "note": "NH₂ auxochrome; free base form."},
         {"lambda_nm": 280, "epsilon": 1400, "solvent": "H₂O pH 11", "note": "B-band enhanced."},
         {"lambda_nm": 203, "epsilon": 7500, "solvent": "H⁺ (protonated)", "note": "Protonation destroys auxochrome effect."},
     ]},
    {"name": "Nitrobenzene", "mw": 123.11, "category": "Nitro Compound",
     "data": [
         {"lambda_nm": 320, "epsilon": 125, "solvent": "hexane", "note": "n→π* weak band."},
         {"lambda_nm": 268, "epsilon": 7800, "solvent": "hexane", "note": "π→π* strong band."},
     ]},
    {"name": "Styrene", "mw": 104.15, "category": "Alkene/Aromatic",
     "data": [
         {"lambda_nm": 248, "epsilon": 14500, "solvent": "hexane", "note": "Vinyl extends conjugation."},
         {"lambda_nm": 282, "epsilon": 450, "solvent": "hexane", "note": "Fine structure shoulder."},
     ]},
    {"name": "Naphthalene", "mw": 128.17, "category": "PAH",
     "data": [
         {"lambda_nm": 275, "epsilon": 5600, "solvent": "ethanol", "note": "α-band (allowed)."},
         {"lambda_nm": 312, "epsilon": 320, "solvent": "ethanol", "note": "p-band (forbidden)."},
         {"lambda_nm": 220, "epsilon": 110000, "solvent": "ethanol", "note": "β-band (very strong)."},
     ]},
    {"name": "Anthracene", "mw": 178.23, "category": "PAH",
     "data": [
         {"lambda_nm": 356, "epsilon": 7900, "solvent": "ethanol", "note": "Visible region onset; fluoresces."},
         {"lambda_nm": 256, "epsilon": 200000, "solvent": "ethanol", "note": "Very intense β-band."},
     ]},
    {"name": "Acetone", "mw": 58.08, "category": "Ketone",
     "data": [
         {"lambda_nm": 279, "epsilon": 15, "solvent": "hexane", "note": "n→π* (forbidden, very weak)."},
         {"lambda_nm": 190, "epsilon": 1000, "solvent": "hexane", "note": "π→π* (VUV)."},
     ]},
    {"name": "Methyl ethyl ketone (MEK)", "mw": 72.11, "category": "Ketone",
     "data": [
         {"lambda_nm": 278, "epsilon": 16, "solvent": "hexane", "note": "n→π* similar to acetone."},
     ]},
    {"name": "Acetaldehyde", "mw": 44.05, "category": "Aldehyde",
     "data": [
         {"lambda_nm": 290, "epsilon": 13, "solvent": "hexane", "note": "n→π* (weak)."},
         {"lambda_nm": 180, "epsilon": 10000, "solvent": "vapor", "note": "π→π*."},
     ]},

    # ═══ Dyes & Chromophores ═══
    {"name": "Methyl Orange", "mw": 327.33, "category": "Dye/Indicator",
     "data": [
         {"lambda_nm": 504, "epsilon": 42000, "solvent": "H₂O", "note": "Acid form (red); azo dye."},
         {"lambda_nm": 464, "epsilon": 25000, "solvent": "H₂O pH>4.4", "note": "Base form (yellow)."},
     ]},
    {"name": "Phenolphthalein", "mw": 318.32, "category": "Indicator",
     "data": [
         {"lambda_nm": 553, "epsilon": 3800, "solvent": "H₂O pH>10", "note": "Basic pink form."},
         {"lambda_nm": "none", "epsilon": 0, "solvent": "H₂O pH<8.2", "note": "Colorless acidic form."},
     ]},
    {"name": "Bromocresol Green", "mw": 698.02, "category": "Indicator",
     "data": [
         {"lambda_nm": 615, "epsilon": 15000, "solvent": "H₂O pH>4.7", "note": "Basic blue form."},
         {"lambda_nm": 444, "epsilon": 7200, "solvent": "H₂O pH<3.8", "note": "Acid yellow form."},
     ]},
    {"name": "Coomassie Brilliant Blue G-250", "mw": 854.03, "category": "Dye",
     "data": [
         {"lambda_nm": 595, "epsilon": 84000, "solvent": "acidic aq", "note": "Bradford protein assay dye-protein complex."},
         {"lambda_nm": 470, "epsilon": 14000, "solvent": "acidic aq", "note": "Free dye (red)."},
     ]},
    {"name": "Fluorescein", "mw": 332.31, "category": "Fluorophore",
     "data": [
         {"lambda_nm": 494, "epsilon": 93000, "solvent": "pH>8 aq", "note": "Excitation max; intense green fluorescence."},
     ]},
    {"name": "Rhodamine B", "mw": 479.02, "category": "Fluorophore/Dye",
     "data": [
         {"lambda_nm": 554, "epsilon": 116000, "solvent": "ethanol", "note": "Extremely intense; fluorescent."},
     ]},

    # ═══ Pharmaceuticals ═══
    {"name": "Aspirin (Acetylsalicylic acid)", "mw": 180.16, "category": "Drug",
     "data": [
         {"lambda_nm": 229, "epsilon": 8500, "solvent": "ethanol", "note": "Ester + aromatic conjugation."},
         {"lambda_nm": 237, "epsilon": 5000, "solvent": "buffer pH 2", "note": "Slightly shifted in acid."},
     ]},
    {"name": "Ibuprofen", "mw": 206.29, "category": "Drug (NSAID)",
     "data": [
         {"lambda_nm": 222, "epsilon": 14500, "solvent": "ethanol/H₂O", "note": "Isobutylphenyl group; aromatic π→π*."},
         {"lambda_nm": 264, "epsilon": 900, "solvent": "ethanol/H₂O", "note": "Weaker B-band."},
     ]},
    {"name": "Paracetamol (Acetaminophen)", "mw": 151.16, "category": "Drug (Analgesic)",
     "data": [
         {"lambda_nm": 245, "epsilon": "~9000", "solvent": "ethanol/0.1N HCl", "note": "Phenol-acetanilide chromophore."},
         {"lambda_nm": 257, "epsilon": "~6000", "solvent": "pH 10 buffer", "note": "Phenolate shift at high pH."},
     ]},
    {"name": "Caffeine", "mw": 194.19, "category": "Drug (Stimulant)",
     "data": [
         {"lambda_nm": 273, "epsilon": 9700, "solvent": "H₂O", "note": "Purine ring system; common HPLC IS."},
     ]},
    {"name": "Doxycycline", "mw": 444.44, "category": "Drug (Antibiotic)",
     "data": [
         {"lambda_nm": 269, "epsilon": 21000, "solvent": "methanol", "note": "Conjugated tetracycline system."},
         {"lambda_nm": 349, "epsilon": 11700, "solvent": "methanol", "note": "Extended chromophore."},
     ]},
    {"name": "Tetracycline", "mw": 444.44, "category": "Drug (Antibiotic)",
     "data": [
         {"lambda_nm": 268, "epsilon": 18000, "solvent": "0.1N HCl", "note": "Acidic form."},
         {"lambda_nm": 355, "epsilon": 12300, "solvent": "0.1N HCl", "note": "Long-wavelength band."},
     ]},
    {"name": "Chloramphenicol", "mw": 323.13, "category": "Drug (Antibiotic)",
     "data": [
         {"lambda_nm": 278, "epsilon": 300, "solvent": "H₂O", "note": "Nitrobenzene-type weak absorption."},
     ]},

    # ═══ Vitamins ═══
    {"name": "Vitamin A (Retinol)", "mw": 286.46, "category": "Vitamin",
     "data": [
         {"lambda_nm": 325, "epsilon": 52000, "solvent": "ethanol", "note": "Polyene chain; highly conjugated."},
         {"lambda_nm": 350, "epsilon": 35000, "solvent": "ethanol", "note": "cis-isomer shifted."},
     ]},
    {"name": "Vitamin B1 (Thiamine)", "mw": 337.27, "category": "Vitamin",
     "data": [
         {"lambda_nm": 246, "epsilon": 12000, "solvent": "pH 7 buffer", "note": "Pyrimidine-pyrazolium system."},
         {"lambda_nm": 235, "epsilon": 11000, "solvent": "0.1N HCl", "note": "Acidic form."},
     ]},
    {"name": "Vitamin B2 (Riboflavin)", "mw": 376.37, "category": "Vitamin",
     "data": [
         {"lambda_nm": 267, "epsilon": 31000, "solvent": "pH 7 buffer", "note": "Isoalloxazine ring."},
         {"lambda_nm": 373, "epsilon": 10000, "solvent": "pH 7 buffer", "note": "Lower energy band."},
         {"lambda_nm": 445, "epsilon": 8200, "solvent": "pH 7 buffer", "note": "Yellow color origin."},
     ]},
    {"name": "Vitamin B12 (Cobalamin)", "mw": 1355.37, "category": "Vitamin",
     "data": [
         {"lambda_nm": 361, "epsilon": 28000, "solvent": "H₂O", "note": "γ-peak; corrin ring."},
         {"lambda_nm": 551, "epsilon": 8300, "solvent": "H₂O", "note": "β-peak; characteristic red color."},
         {"lambda_nm": 528, "epsilon": 7000, "solvent": "H₂O", "note": "α-peak."},
     ]},
    {"name": "Vitamin C (Ascorbic acid)", "mw": 176.12, "category": "Vitamin",
     "data": [
         {"lambda_nm": 265, "epsilon": 14500, "solvent": "pH 3 buffer", "note": "Enediol lactone structure."},
         {"lambda_nm": 244, "epsilon": 11000, "solvent": "0.01N H₂SO₄", "note": "Acidic form."},
     ]},

    # ═══ Pesticides / Environmental ═══
    {"name": "DDT", "mw": 354.49, "category": "Pesticide (OC)",
     "data": [
         {"lambda_nm": 236, "epsilon": 12400, "solvent": "hexane", "note": "Trichloroethyl-aromatic."},
     ]},
    {"name": "Atrazine", "mw": 215.68, "category": "Pesticide (Herbicide)",
     "data": [
         {"lambda_nm": 220, "epsilon": 34000, "solvent": "methanol", "note": "Triazine ring + chloro substituents."},
         {"lambda_nm": 265, "epsilon": 4000, "solvent": "methanol", "note": "Weaker n→π* band."},
     ]},
    {"name": "Carbofuran", "mw": 222.26, "category": "Pesticide (Carbamate)",
     "data": [
         {"lambda_nm": 278, "epsilon": 2500, "solvent": "methanol", "note": "Carbamate-furan system."},
     ]},

    # ═══ Metal Complexes / Coordination Compounds ═══
    {"name": "[Fe(phen)₃]²⁺ (Ferroin)", "mw": 696.51, "category": "Metal Complex",
     "data": [
         {"lambda_nm": 510, "epsilon": 11100, "solvent": "H₂O", "note": "MLCT band; intense orange-red color. Redox indicator."},
     ]},
    {"name": "[Cu(NH₃)₄]²⁺", "mw": 227.74, "category": "Metal Complex",
     "data": [
         {"lambda_nm": 600, "epsilon": 120, "solvent": "aq excess NH₃", "note": "d-d transition; deep blue color (weak but visible)."},
     ]},
    {"name": "KMnO₄", "mw": 158.03, "category": "Inorganic Salt",
     "data": [
         {"lambda_nm": 525, "epsilon": 2230, "solvent": "H₂O", "note": "LMCT; intense purple color."},
         {"lambda_nm": 545, "epsilon": 2200, "solvent": "H₂O", "note": "Shoulder."},
     ]},
    {"name": "K₂Cr₂O₇", "mw": 294.18, "category": "Inorganic Salt",
     "data": [
         {"lambda_nm": 350, "epsilon": 1716, "solvent": "0.05M H₂SO₄", "note": "LMCT; orange color. Standard for spectrophotometer calibration."},
         {"lambda_nm": 440, "epsilon": 15, "solvent": "0.05M H₂SO₄", "note": "Weak d-d tail."},
         {"lambda_nm": 257, "epsilon": 14400, "solvent": "0.05M H₂SO₄", "note": "Charge transfer (strong)."},
     ]},
]


@ChemMCPManager.register_tool
class MolarAbsorptivityLookup(BaseTool):
    """
    摩尔吸光系数查询工具：从内置数据库中查找化合物在特定波长下的 ε 值，
    支持按化合物名称、类别、波长范围检索。
    """
    __version__ = "0.1.0"
    name = "MolarAbsorptivityLookup"
    func_name = "lookup_molar_absorptivity"
    description = "Look up molar absorptivity (ε, M⁻¹cm⁻¹) values from a database of common compounds across various wavelengths and solvents."
    implementation_description = "Built-in database of ~45 compounds including amino acids, nucleotides, pharmaceuticals, dyes, vitamins, pesticides, and metal complexes. Supports search by name, category, wavelength range, or minimum ε threshold."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Molar Absorptivity", "UV-Vis", "Database", "Spectroscopy", "Beer-Lambert", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("compound_name", "str", "", "Compound name to look up (exact or fuzzy match)."),
        ("category", "str", "", "Filter by category (e.g., 'Amino Acid', 'Drug', 'Dye'). Leave empty for all."),
        ("wavelength_nm", "float", "0", "Find compounds with data near this wavelength (±10nm tolerance)."),
        ("min_epsilon", "float", "0", "Minimum molar absorptivity value filter."),
        ("max_epsilon", "float", "0", "Maximum molar absorptivity filter (0 = no limit)."),
        ("solvent", "str", "", "Filter by solvent (e.g., 'water', 'ethanol', 'hexane')."),
        ("search_mode", "str", "name", "Search mode: 'name', 'category', 'wavelength', 'browse_all'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "E.g., 'tryptophan' or 'wavelength 280' or 'category Drug'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with matching entries, ε values, wavelengths, solvents, and usage notes."),
    ]

    examples = [
        {
            "code_input": {
                "compound_name": "tryptophan",
                "category": "",
                "wavelength_nm": 0,
                "min_epsilon": 0,
                "max_epsilon": 0,
                "solvent": "",
                "search_mode": "name",
            },
            "text_input": {
                "input_params": "tryptophan",
            },
            "output": {
                "result": {
                    "mode": "lookup",
                    "compound": "Tryptophan",
                    "entries_found": 1,
                    "note": "Database match result.",
                }
            }
        },
        {
            "code_input": {
                "compound_name": "",
                "category": "",
                "wavelength_nm": 280,
                "min_epsilon": 1000,
                "max_epsilon": 0,
                "solvent": "",
                "search_mode": "wavelength",
            },
            "text_input": {
                "input_params": "wavelength 280",
            },
            "output": {
                "result": {
                    "mode": "wavelength_search",
                    "wavelength_nm": 280,
                    "tolerance_nm": 10,
                    "note": "Compounds absorbing near 280nm.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _search_by_name(self, name: str) -> List[dict]:
        name_lower = name.lower().strip()
        results = []
        for entry in EPSILON_DB:
            if name_lower in entry["name"].lower() or entry["name"].lower() in name_lower:
                results.append(entry)
        return results

    def _search_by_category(self, category: str) -> List[dict]:
        cat_lower = category.lower().strip()
        results = []
        for entry in EPSILON_DB:
            if cat_lower in entry["category"].lower():
                results.append(entry)
        return results

    def _search_by_wavelength(self, wl: float, min_eps: float = 0,
                               max_eps: float = 0, solvent: str = "") -> List[dict]:
        tolerance = 10  # nm
        results = []
        for entry in EPSILON_DB:
            for dp in entry.get("data", []):
                lam = dp["lambda_nm"]
                eps = dp["epsilon"]
                sol = dp.get("solvent", "")
                if abs(lam - wl) <= tolerance:
                    if min_eps > 0 and eps < min_eps:
                        continue
                    if max_eps > 0 and eps > max_eps:
                        continue
                    if solvent and solvent.lower() not in sol.lower():
                        continue
                    matched = {**entry, "matched_data_point": dp}
                    if matched not in results:
                        results.append(matched)
        return results

    def _run_base(self, compound_name: str = "", category: str = "",
                  wavelength_nm: float = 0.0, min_epsilon: float = 0.0,
                  max_epsilon: float = 0.0, solvent: str = "",
                  search_mode: str = "name") -> dict:

        mode = search_mode.lower().strip()

        if mode == "wavelength" and wavelength_nm > 0:
            results = self._search_by_wavelength(
                wavelength_nm, min_epsilon, max_epsilon, solvent)
            summary = f"Wavelength search near {wavelength_nm}±10nm"
        elif mode == "category" and category:
            results = self._search_by_category(category)
            summary = f"Category search: '{category}'"
        elif mode == "browse_all":
            results = EPSILON_DB
            summary = "Browse all entries"
        elif compound_name:
            results = self._search_by_name(compound_name)
            summary = f"Name search: '{compound_name}'"
        else:
            raise ChemMCPError(
                "Provide either compound_name, category, wavelength_nm (>0), "
                "or search_mode='browse_all'."
            )

        if not results:
            return {"result": {
                "mode": "lookup",
                "search_summary": summary,
                "entries_found": 0,
                "message": "No matches found. Try different search terms.",
                "available_categories": sorted(set(e["category"] for e in EPSILON_DB)),
                "total_database_entries": len(EPSILON_DB),
            }}

        # Format results
        formatted = []
        for r in results:
            entry = {
                "name": r["name"],
                "mw": r["mw"],
                "category": r["category"],
            }
            if "matched_data_point" in r:
                dp = r["matched_data_point"]
                entry["matched_wavelength_nm"] = dp["lambda_nm"]
                entry["molar_absorptivity"] = dp["epsilon"]
                entry["solvent"] = dp.get("solvent", "")
                entry["note"] = dp.get("note", "")
            else:
                entry["all_wavelength_data"] = r.get("data", [])
            formatted.append(entry)

        return {"result": {
            "mode": "lookup",
            "search_summary": summary,
            "entries_found": len(formatted),
            "matches": formatted,
            "database_info": {
                "total_compounds": len(EPSILON_DB),
                "categories_available": sorted(set(e["category"] for e in EPSILON_DB)),
                "wavelength_range_nm": [190, 615],
                "epsilon_range": [0, 200000],
            },
            "usage_tip": (
                "Use ε values in Beer-Lambert law: A = ε·b·c. "
                "For quantitative analysis, choose λ where ε is large and the signal is linear."
            ),
        }}

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            first = parts[0].lower() if parts else ""
            if first == "wavelength" and len(parts) > 1:
                wl = float(parts[1])
                return self._run_base(wavelength_nm=wl, search_mode="wavelength")
            elif first == "category" and len(parts) > 1:
                cat = " ".join(parts[1:])
                return self._run_base(category=cat, search_mode="category")
            elif first == "all" or first == "browse":
                return self._run_base(search_mode="browse_all")
            else:
                return self._run_base(compound_name=input_params.strip(), search_mode="name")
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input '{input_params}': {e}")
