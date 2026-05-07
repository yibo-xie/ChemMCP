"""
基组处理工具 (Basis Set Handler) — MCP #461
处理常见量子化学基组，支持 STO/GTO 转换、基组信息查询与比较。
"""
import logging
import math
from typing import Optional, List, Dict, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BasisSetHandler(BaseTool):
    """
    基组处理工具。支持 STO-nG, 3-21G, 6-31G*, 6-311+G**, cc-pVDZ 等常见基组的
    信息查询、STO↔GTO 展开系数转换、基组大小/收缩度分析。
    """
    __version__ = "0.1.0"
    name = "BasisSetHandler"
    func_name = "basis_set_handler"
    description = "Handle quantum chemistry basis sets: query info, convert STO to GTO expansion coefficients, compare basis sets, analyze contraction/size."
    implementation_description = "Implements a basis set database for common quantum chemistry basis sets (STO-nG, Pople-style, Dunning correlation-consistent). Provides STO→GTO expansion coefficients, contraction schemes, and basis set comparison metrics."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Chemistry", "Basis Set", "STO", "GTO", "Basis Set Conversion", "Electronic Structure"]
    required_envs = []

    code_input_sig = [
        ("basis_set_name", "str", "N/A", "Basis set name: 'STO-3G', 'STO-6G', '3-21G', '6-31G', '6-31G*', '6-311G**', '6-311++G**', 'cc-pVDZ', 'cc-pVTZ', 'aug-cc-pVDZ', 'def2-SVP', 'def2-TZVP'."),
        ("operation", "str", "'info'", "Operation: 'info' (query details), 'sto_gto_convert' (get expansion), 'compare' (compare two sets), 'list' (all available), 'contraction' (scheme analysis)."),
        ("element", "str", "'H'", "Element symbol for element-specific data (case-insensitive)."),
        ("n_primitives", "int", "3", "Number of primitive GTOs for STO-nG conversion (n in STO-nG)."),
        ("compare_with", "str", "None", "Second basis set name for comparison (used with operation='compare')."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Space-separated: basis_set_name operation [element] [n_primitives] [compare_with]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing basis set information, expansion coefficients, comparison results, or list of available sets."),
    ]

    examples = [
        {
            "code_input": {
                "basis_set_name": "STO-3G",
                "operation": "info",
                "element": "C",
            },
            "text_input": {
                "input_str": "STO-3G info C",
            },
            "output": {
                "result": {
                    "basis_set": "STO-3G",
                    "element": "C",
                    "type": "Minimal Basis",
                    "n_primitives_per_orbital": 3,
                    "description": "Slater-Type Orbital approximated by 3 Gaussian functions",
                }
            }
        },
        {
            "code_input": {
                "basis_set_name": "STO-3G",
                "operation": "sto_gto_convert",
                "element": "H",
                "n_primitives": 3,
            },
            "text_input": {
                "input_str": "STO-3G sto_gto_convert H 3",
            },
            "output": {
                "result": {
                    "expansion_coefficients": [...],
                    "exponents": [...],
                }
            }
        },
    ]

    # ── Basis Set Database ──────────────────────────────────────────
    # Standard STO-nG exponents and coefficients (Hehre-Stewart-Pople)
    # Reference: Hehre, W.J., Stewart, R.F., Pople, J.A. (1969/1972)
    _STO_NG_DATA = {
        # element: {orbital: [(d_i, α_i) for each primitive], ...}
        "H": {
            "1s": {
                2: [  # STO-2G
                    (0.430129, 0.130976),
                    (0.678914, 0.478319),
                    (1.0, 0.0),  # placeholder normalization handled separately
                ],
                3: [  # STO-3G
                    (0.154329, 0.109818),
                    (0.535328, 0.405771),
                    (0.444635, 2.227660),
                ],
                4: [  # STO-4G
                    (0.019685, 0.049282),
                    (0.137965, 0.234827),
                    (0.478319, 0.834525),
                    (0.363812, 3.025230),
                ],
                5: [  # STO-5G
                    (0.009163, 0.027707),
                    (0.049451, 0.138467),
                    (0.168472, 0.424057),
                    (0.370747, 0.785842),
                    (0.294602, 2.615730),
                ],
                6: [  # STO-6G
                    (0.004537, 0.016653),
                    (0.029644, 0.094069),
                    (0.117885, 0.310385),
                    (0.320283, 0.656480),
                    (0.407858, 1.247790),
                    (0.279794, 3.625680),
                ],
            },
            "zeta_1s": 1.0,  # orbital exponent ζ for H 1s
        },
        "C": {
            "1s": {
                3: [
                    (0.154329, 0.072586),
                    (0.535328, 0.291208),
                    (0.444635, 1.242567),
                ],
                4: [
                    (0.019685, 0.032592),
                    (0.137965, 0.166159),
                    (0.478319, 0.514008),
                    (0.363812, 1.847951),
                ],
                6: [
                    (0.004537, 0.010875),
                    (0.029644, 0.061455),
                    (0.117885, 0.233202),
                    (0.320283, 0.532841),
                    (0.407858, 1.037350),
                    (0.279794, 3.119150),
                ],
            },
            "2s": {
                3: [
                    (-0.099967, 0.072586),
                    (0.399513, 0.291208),
                    (0.700115, 1.242567),
                ],
                4: [
                    (-0.016674, 0.032592),
                    (0.097494, 0.166159),
                    (0.577201, 0.514008),
                    (0.341981, 1.847951),
                ],
                6: [
                    (-0.003519, 0.010875),
                    (-0.023114, 0.061455),
                    (-0.134003, 0.233202),
                    (0.468857, 0.532841),
                    (0.601685, 1.037350),
                    (0.245179, 3.119150),
                ],
            },
            "2p": {
                3: [
                    (0.155916, 0.072586),
                    (0.607684, 0.291208),
                    (0.391957, 1.242567),
                ],
                4: [
                    (0.022766, 0.032592),
                    (0.151268, 0.166159),
                    (0.503411, 0.514008),
                    (0.322255, 1.847951),
                ],
                6: [
                    (0.004912, 0.010875),
                    (0.033867, 0.061455),
                    (0.173972, 0.233202),
                    (0.518176, 0.532841),
                    (0.455425, 1.037350),
                    (0.167655, 3.119150),
                ],
            },
            "zeta_1s": 5.67,   # C 1s orbital exponent
            "zeta_2sp": 1.72,  # C 2sp orbital exponent (same for 2s and 2p in STO-nG)
        },
        "N": {
            "1s": {3: [(0.154329, 0.098766), (0.535328, 0.393536), (0.444635, 1.643667)]},
            "2s": {3: [(-0.099967, 0.098766), (0.399513, 0.393536), (0.700115, 1.643667)]},
            "2p": {3: [(0.155916, 0.098766), (0.607684, 0.393536), (0.391957, 1.643667)]},
            "zeta_1s": 6.67, "zeta_2sp": 1.95,
        },
        "O": {
            "1s": {3: [(0.154329, 0.127817), (0.535328, 0.467871), (0.444635, 2.108810)]},
            "2s": {3: [(-0.099967, 0.127817), (0.399513, 0.467871), (0.700115, 2.108810)]},
            "2p": {3: [(0.155916, 0.127817), (0.607684, 0.467871), (0.391957, 2.108810)]},
            "zeta_1s": 7.66, "zeta_2sp": 2.25,
        },
        "F": {
            "1s": {3: [(0.154329, 0.160856), (0.535328, 0.545813), (0.444635, 2.618710)]},
            "2s": {3: [(-0.099967, 0.160856), (0.399513, 0.545813), (0.700115, 2.618710)]},
            "2p": {3: [(0.155916, 0.160856), (0.607684, 0.545813), (0.391957, 2.618710)]},
            "zeta_1s": 8.66, "zeta_2sp": 2.55,
        },
        "He": {
            "1s": {3: [(0.154329, 0.231266), (0.535328, 0.836218), (0.444635, 3.850620)]},
            "zeta_1s": 1.6875,
        },
        "Li": {
            "1s": {3: [(0.154329, 0.162205), (0.535328, 0.606796), (0.444635, 3.206570)]},
            "2s": {3: [(-0.099967, 0.162205), (0.399513, 0.606796), (0.700115, 3.206570)]},
            "2p": {3: [(0.155916, 0.162205), (0.607684, 0.606796), (0.391957, 3.206570)]},
            "zeta_1s": 2.69, "zeta_2sp": 0.65,
        },
        "Be": {
            "1s": {3: [(0.154329, 0.212024), (0.535328, 0.770509), (0.444635, 2.945260)]},
            "2s": {3: [(-0.099967, 0.212024), (0.399513, 0.770509), (0.700115, 2.945260)]},
            "2p": {3: [(0.155916, 0.212024), (0.607684, 0.770509), (0.391957, 2.945260)]},
            "zeta_1s": 3.69, "zeta_2sp": 0.975,
        },
        "B": {
            "1s": {3: [(0.154329, 0.267058), (0.535328, 0.935511), (0.444635, 2.741830)]},
            "2s": {3: [(-0.099967, 0.267058), (0.399513, 0.935511), (0.700115, 2.741830)]},
            "2p": {3: [(0.155916, 0.267058), (0.607684, 0.935511), (0.391957, 2.741830)]},
            "zeta_1s": 4.68, "zeta_2sp": 1.30,
        },
    }

    # ── Basis Set Metadata Database ────────────────────────────────
    _BASIS_INFO = {
        "STO-3G": {
            "type": "Minimal Basis",
            "category": "STO-nG",
            "n_primitives": 3,
            "description": "Each STO fitted by 3 GTOs. Minimal basis (one function per valence orbital).",
            "accuracy": "Qualitative",
            "typical_error_eV": "~1-2 eV",
            "functions_per_atom_H": 1,
            "functions_per_atom_C": 5,  # 1s + 2s + 2px + 2py + 2pz
            "polarization": False,
            "diffuse": False,
            "year": 1969,
            "reference": "Hehre, Stewart & Pople, JCP 1969",
        },
        "STO-6G": {
            "type": "Minimal Basis",
            "category": "STO-nG",
            "n_primitives": 6,
            "description": "Each STO fitted by 6 GTOs. Better accuracy than STO-3G.",
            "accuracy": "Semi-quantitative",
            "typical_error_eV": "~0.5-1 eV",
            "functions_per_atom_H": 1,
            "functions_per_atom_C": 5,
            "polarization": False,
            "diffuse": False,
            "year": 1969,
        },
        "3-21G": {
            "type": "Split Valence",
            "category": "Pople",
            "n_primitives_split": "3 inner / 2 outer",
            "description": "Valence orbitals split into 2 contracted functions (3 primitives → 1 inner + 2 outer). Core: 3 primitives.",
            "accuracy": "Semi-quantitative",
            "typical_error_eV": "~0.3-0.8 eV",
            "functions_per_atom_H": 2,
            "functions_per_atom_C": 9,
            "polarization": False,
            "diffuse": False,
            "year": 1980,
            "reference": "Binkley, Pople & Melius, JCP 1980",
        },
        "6-31G": {
            "type": "Split Valence",
            "category": "Pople",
            "n_primitives_split": "6 inner / 3 outer (valence)",
            "description": "Valence split into 2 contractions (6+3 primitives). Better than 3-21G.",
            "accuracy": "Good",
            "typical_error_eV": "~0.2-0.5 eV",
            "functions_per_atom_H": 2,
            "functions_per_atom_C": 9,
            "polarization": False,
            "diffuse": False,
            "year": 1984,
            "reference": "Ditchfield, Hehre & Pople, JCP 1971; Hariharan & Pople, Theo Chim Acta 1973",
        },
        "6-31G*": {
            "type": "Split Valence + Polarization",
            "category": "Pople",
            "n_primitives_split": "6 inner / 3 outer (valence) + d on heavy atoms",
            "description": "6-31G plus d-type polarization functions on heavy atoms (Z>2).",
            "accuracy": "Very Good",
            "typical_error_eV": "~0.1-0.3 eV",
            "functions_per_atom_H": 2,
            "functions_per_atom_C": 15,  # 9 + 6 d-functions
            "polarization": True,
            "polarization_functions": "d (5 or 6 components) on heavy atoms",
            "diffuse": False,
            "year": 1986,
            "reference": "Pietro et al., JCP 1982; Harharan & Pople, Theo Chim Acta 1973",
        },
        "6-31G**": {
            "type": "Split Valence + Full Polarization",
            "category": "Pople",
            "description": "6-31G* plus p-type polarization on H as well.",
            "accuracy": "Very Good",
            "functions_per_atom_H": 4,  # 1s + 3 p
            "functions_per_atom_C": 18,  # 15 + 3 p on H equivalent
            "polarization": True,
            "polarization_functions": "d on heavy atoms, p on H",
            "diffuse": False,
            "year": 1986,
        },
        "6-311G**": {
            "type": "Triple Split Valence + Full Polarization",
            "category": "Pople",
            "description": "Valence split into 3 contractions (6+3+1). Plus polarization.",
            "accuracy": "Excellent",
            "functions_per_atom_H": 5,
            "functions_per_atom_C": 24,
            "polarization": True,
            "diffuse": False,
            "year": 1987,
            "reference": "Krishnan et al., JCP 1980",
        },
        "6-311++G**": {
            "type": "Triple Split + Full Polarization + Diffuse",
            "category": "Pople",
            "description": "6-311G** plus diffuse s functions on all atoms.",
            "accuracy": "Excellent (anions, weak interactions)",
            "functions_per_atom_H": 6,
            "functions_per_atom_C": 26,
            "polarization": True,
            "diffuse": True,
            "year": 1989,
        },
        "cc-pVDZ": {
            "type": "Correlation Consistent Double Zeta",
            "category": "Dunning",
            "description": "Dunning's cc-pVDZ: systematically convergent to CBS limit. Includes polarization.",
            "accuracy": "Excellent for correlated methods",
            "functions_per_atom_H": 2,
            "functions_per_atom_C": 14,
            "polarization": True,
            "diffuse": False,
            "year": 1989,
            "reference": "Dunning, JCP 1989",
        },
        "cc-pVTZ": {
            "type": "Correlation Consistent Triple Zeta",
            "category": "Dunning",
            "description": "Triple-zeta version of cc-pV series. Higher angular momentum functions.",
            "accuracy": "Near-CBS quality for many properties",
            "functions_per_atom_H": 5,
            "functions_per_atom_C": 30,
            "polarization": True,
            "diffuse": False,
            "year": 1989,
        },
        "aug-cc-pVDZ": {
            "type": "Augmented Correlation Consistent Double Zeta",
            "category": "Dunning",
            "description": "cc-pVDZ plus diffuse functions on all atoms. Excellent for anions/van der Waals.",
            "accuracy": "Excellent (anions, excited states, weak interactions)",
            "functions_per_atom_H": 5,
            "functions_per_atom_C": 22,
            "polarization": True,
            "diffuse": True,
            "year": 1995,
            "reference": "Kendall, Harrison & Dunning, JCP 1992",
        },
        "def2-SVP": {
            "type": "Split Valence Polarized",
            "category": "Ahlrichs",
            "description": "Ahrlighs def2 split-valence polarized. General-purpose efficiency.",
            "accuracy": "Good",
            "functions_per_atom_H": 4,
            "functions_per_atom_C": 17,
            "polarization": True,
            "diffuse": False,
            "year": 2004,
            "reference": "Weigend & Ahlrichs, Phys Chem Chem Phys 2005",
        },
        "def2-TZVP": {
            "type": "Triple Zeta Valence Polarized",
            "category": "Ahlrichs",
            "description": "def2 triple-zeta valence polarized. Higher accuracy than SVP.",
            "accuracy": "Very Good to Excellent",
            "functions_per_atom_H": 5,
            "functions_per_atom_C": 30,
            "polarization": True,
            "diffuse": False,
            "year": 2004,
        },
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.Hartree_to_eV = 27.211386245988

    def _run_base(self, basis_set_name: str, operation: str = "info",
                  element: str = "H", n_primitives: int = 3,
                  compare_with: str = None) -> dict:
        """Core logic."""
        op = operation.lower().strip()
        bs = basis_set_name.strip()

        # Special handling: if basis_set_name looks like an operation command
        _command_like = {"list", "compare", "sto_gto_convert", "convert", "contraction"}
        if bs.lower().strip() in _command_like:
            op = bs.lower().strip()
            bs = "STO-3G"  # placeholder for operations that don't need a specific set
        else:
            op = operation.lower().strip()

        if op == "list":
            return self._list_basis_sets()
        elif op == "info":
            return self._bs_info(bs, element)
        elif op in ("sto_gto_convert", "convert"):
            return self._sto_gto_convert(element, n_primitives)
        elif op == "compare":
            return self._compare_basis_sets(bs, compare_with, element)
        elif op == "contraction":
            return self._contraction_scheme(bs, element)
        else:
            raise ChemMCPError(
                f"Unknown operation '{operation}'. "
                f"Use: info, sto_gto_convert, compare, list, contraction."
            )

    # ── List All Available Basis Sets ──────────────────────────────
    def _list_basis_sets(self) -> dict:
        return {"result": {
            "available_basis_sets": list(self._BASIS_INFO.keys()),
            "total_count": len(self._BASIS_INFO),
            "categories": sorted(set(v["category"] for v in self._BASIS_INFO.values())),
            "supported_elements": list(self._STO_NG_DATA.keys()),
            "note": "For STO-nG expansion coefficients, use sto_gto_convert operation.",
        }}

    # ── Basis Set Info Query ───────────────────────────────────────
    def _bs_info(self, bs: str, el: str) -> dict:
        if bs not in self._BASIS_INFO:
            available = ", ".join(sorted(self._BASIS_INFO.keys()))
            raise ChemMCPError(
                f"Unknown basis set '{bs}'. Available: {available}"
            )
        info = dict(self._BASIS_INFO[bs])
        info["basis_set"] = bs  # include the queried basis set name
        el_upper = el.upper().strip()

        # Add element-specific data
        if el_upper in self._STO_NG_DATA:
            elem_data = self._STO_NG_DATA[el_upper]
            zeta_keys = [k for k in elem_data if k.startswith("zeta_")]
            if zeta_keys:
                info["element_data"] = {
                    "element": el_upper,
                    "orbital_exponents": {k: v for k, v in elem_data.items() if k.startswith("zeta_")},
                    "available_orbitals": [k for k in elem_data.keys() if not k.startswith("zeta")],
                    "available_n_values": sorted(set(
                        n for orb in elem_data.values() if isinstance(orb, dict)
                        for n in orb.keys()
                    )),
                }

        return {"result": info}

    # ── STO → GTO Expansion Coefficients ───────────────────────────
    def _sto_gto_convert(self, element: str, n: int) -> dict:
        el = element.upper().strip()
        if el not in self._STO_NG_DATA:
            raise ChemMCPError(
                f"Element '{el}' not found. Available: {list(self._STO_NG_DATA.keys())}"
            )
        elem_data = self._STO_NG_DATA[el]

        expansions = {}
        for orbital, n_dict in elem_data.items():
            if orbital.startswith("zeta_"):
                continue
            if n not in n_dict:
                continue
            coeffs = n_dict[n]
            # Normalize: Σ d_i = 1 (for properly normalized STO-nG)
            total_d = sum(c[0] for c in coeffs)
            normalized = [(d / total_d, alpha) for d, alpha in coeffs]

            # Compute the actual GTO representation
            # φ_STO(r) ≈ Σ_i d_i · g_i(r) where g_i(r) = N(α_i) exp(-α_i r²)
            expansions[orbital] = {
                "n_primitives": n,
                "raw_coefficients_exponents": [{"d_i": c[0], "alpha_i": c[1]} for c in coeffs],
                "normalized_coefficients": [{"d_i_norm": d, "alpha_i": a} for d, a in normalized],
                "sum_check_d": round(total_d, 6),
                "formula": f"φ_{orbital}^{{STO}}(r) ≈ Σ_{{i=1}}^{{n}} d_i · (2α_i/π)^{{3/4}} exp(-α_i r²)",
                "fitting_quality": "Least-squares fit to Slater-type orbital",
            }

        zeta_info = {k: v for k, v in elem_data.items() if k.startswith("zeta_")}
        return {"result": {
            "element": el,
            "n_gaussians": n,
            "orbital_expansions": expansions,
            "orbital_exponents_zeta": zeta_info,
            "note": f"STO-{n}G: Each Slater-type orbital fitted by {n} primitive Gaussians. "
                   f"Coefficients from least-squares minimization of ∫[φ_STO - φ_GTO]² dr.",
        }}

    # ── Compare Two Basis Sets ─────────────────────────────────────
    def _compare_basis_sets(self, bs1: str, bs2: str, el: str) -> dict:
        if bs1 not in self._BASIS_INFO or bs2 not in self._BASIS_INFO:
            raise ChemMCPError(f"Both basis sets must be valid. Available: {list(self._BASIS_INFO.keys())}")

        info1 = self._BASIS_INFO[bs1]
        info2 = self._BASIS_INFO[bs2]

        comparison = {
            "basis_set_A": bs1,
            "basis_set_B": bs2,
            "comparison_element": el.upper(),
            "size_comparison": {
                "A_functions": info1.get("functions_per_atom_" + el.upper(), "N/A"),
                "B_functions": info2.get("functions_per atom_" + el.upper(), "N/A"),
            },
            "features": {
                "A_polarization": info1.get("polarization", False),
                "B_polarization": info2.get("polarization", False),
                "A_diffuse": info1.get("diffuse", False),
                "B_diffuse": info2.get("diffuse", False),
            },
            "accuracy": {
                "A_accuracy": info1.get("accuracy", "Unknown"),
                "B_accuracy": info2.get("accuracy", "Unknown"),
                "A_typical_error": info1.get("typical_error_eV", "Unknown"),
                "B_typical_error": info2.get("typical_error_eV", "Unknown"),
            },
            "recommendation": self._generate_recommendation(info1, info2),
        }
        return {"result": comparison}

    # ── Contraction Scheme Analysis ────────────────────────────────
    def _contraction_scheme(self, bs: str, el: str) -> dict:
        if bs not in self._BASIS_INFO:
            raise ChemMCPError(f"Unknown basis set '{bs}'.")

        info = self._BASIS_INFO[bs]
        el_upper = el.upper().strip()

        scheme = {
            "basis_set": bs,
            "element": el_upper,
            "type": info["type"],
            "category": info["category"],
        }

        # Build detailed contraction scheme based on basis type
        if bs.startswith("STO-"):
            n = int(bs.split("-")[1].replace("G", ""))
            scheme.update({
                "contraction_pattern": f"STO-{n}G: {n} primitive GTOs → 1 contracted function per AO",
                "example_contraction": self._get_sto_ng_contraction(el_upper, n),
            })
        elif bs.startswith("6-31"):
            scheme.update({
                "contraction_pattern": "Core: 6 primitives → 1 function; Valence: 3 primitives split into 2 contractions",
                "core_contraction": "(6) → 1 function",
                "valence_contraction": "(3) → 2 functions (inner + outer)",
            })
            if "*" in bs:
                scheme["polarization"] = "Uncontracted d functions added on heavy atoms"
                if "**" in bs:
                    scheme["polarization_h"] = "Uncontracted p functions added on H"
            if "++" in bs:
                scheme["diffuse"] = "Diffuse s (and p) functions added"
        elif bs.startswith("cc-p"):
            scheme.update({
                "contraction_pattern": "Correlation-consistent hierarchical contraction",
                "principle": "Functions grouped by angular momentum, optimized for systematic convergence to CBS limit",
            })
            if "aug-" in bs:
                scheme["augmentation"] = "Additional diffuse shell of each angular momentum"

        return {"result": scheme}

    # ── Helpers ────────────────────────────────────────────────────
    @staticmethod
    def _get_sto_ng_contraction(el: str, n: int) -> str:
        if el == "H":
            return f"H 1s: ({n}) → [1]  (e.g., {n} primitives contracted to one 1s function)"
        elif el == "C":
            return f"C 1s: ({n}) → [1]; C 2s: ({n}) → [1]; C 2px/py/pz: ({n}) → [1]"
        else:
            return f"{el}: ({n}) → [1] per atomic orbital"

    @staticmethod
    def _generate_recommendation(info1: dict, info2: dict) -> str:
        acc_order = {"Qualitative": 0, "Semi-quantitative": 1, "Good": 2,
                     "Very Good": 3, "Excellent": 4, "Near-CBS quality": 5}
        a1 = acc_order.get(info1.get("accuracy", ""), -1)
        a2 = acc_order.get(info2.get("accuracy", ""), -1)

        if a1 > a2:
            return f"{info1.get('type','')} generally provides better accuracy than {info2.get('type','')}. "
        elif a2 > a1:
            return f"{info2.get('type','')} generally provides better accuracy than {info1.get('type','')}. "
        else:
            return "Both basis sets have similar accuracy levels. Choice depends on specific application needs."

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.strip().split()
            bs = parts[0]
            op = parts[1] if len(parts) > 1 else "info"
            el = parts[2] if len(parts) > 2 else "H"
            n = int(parts[3]) if len(parts) > 3 else 3
            cmp = parts[4] if len(parts) > 4 else None
            return self._run_base(bs, op, el, n, cmp)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
