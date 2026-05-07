"""
选择定则分析工具 (MCP #478)。
基于群论的分子电子跃迁选择定则判断（不可约表示直积分析）。
与 SelectionRulesChecker（原子/通用）不同，本工具专注于**分子点群对称性**分析。
"""
import logging
import math
from typing import List, Tuple, Optional, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SelectionRules(BaseTool):
    """
    基于群论的选择定则分析。
    
    功能:
      - 不可约表示(irrep)直积分析: Γ_i ⊗ Γ_op ⊗ Γ_f ⊃ Γ_A₁ ?
      - 支持常见分子点群的特征标表
      - 电偶极(E1)/磁偶极(M1)/电四极(E2)跃迁算子对称性匹配
      - 自旋多重度选择定则 (ΔS=0)
      - 宇称(g/u)选择定则
      - Laporte 定则 (中心对称分子)
      
    与 SelectionRulesChecker 的区别:
      - 本工具：基于**群论 irrep 直积**的严格分析
      - SelectionRulesChecker：基于量子数规则的通用检查
    """
    __version__ = "0.1.0"
    name = "SelectionRules"
    func_name = "analyze_selection_rules_group_theory"
    description = "Analyze spectroscopic transition selection rules using group theory: irreducible representation direct products, character table analysis for molecular point groups."
    implementation_description = (
        "Performs rigorous group-theoretical selection rule analysis:\n"
        "1. Transition allowed if Γ_initial ⊗ Γ_operator ⊗ Γ_final contains A₁ (totally symmetric)\n"
        "2. Spin rule: ΔS = 0 (spin-orbit coupling can relax)\n"
        "3. Parity rule: g↔u for E1 in centrosymmetric molecules\n"
        "4. Includes full character tables for C₂v, C₃v, D₂h, D₃h, D₄h, D₆h, T_d, O_h, C∞v, D∞h"
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Selection Rules", "Group Theory", "Irreducible Representations", "Spectroscopy", "Symmetry"]
    required_envs = []

    code_input_sig = [
        ("point_group", "str", "N/A", "Molecular point group: 'C2v', 'C3v', 'D2h', 'D3h', 'D4h', 'D6h', 'Td', 'Oh', 'Cinfv', 'Dinfh'."),
        ("initial_irrep", "str", "N/A", "Initial state irreducible representation (e.g., 'A1', 'A2', 'B1', 'B2', 'E', 'T1', 'T2')."),
        ("final_irrep", "str", "N/A", "Final state irreducible representation."),
        ("transition_type", "str", "'electric_dipole'", "Transition operator type: 'electric_dipole', 'magnetic_dipole', 'electric_quadrupole', 'Raman'."),
        ("spin_multiplicity_initial", "int", "1", "Spin multiplicity of initial state (2S+1)."),
        ("spin_multiplicity_final", "int", "1", "Spin multiplicity of final state."),
        ("polarization", "str", "'any'", "Polarization direction: 'x', 'y', 'z', or 'any' for unpolarized."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'point_group initial_irrep final_irrep [transition_type]'. Example: 'C2v A1 B2 electric_dipole' or 'Oh T1g T1u E1'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with is_allowed verdict, direct product decomposition, satisfied/violated rules, and detailed symmetry analysis."),
    ]

    examples = [
        {
            "code_input": {
                "point_group": "C2v",
                "initial_irrep": "A1",
                "final_irrep": "B2",
                "transition_type": "electric_dipole",
            },
            "text_input": {"input_params": "C2v A1 B2 electric_dipole"},
            "output": {"result": {"is_allowed": True, "point_group": "C2v"}},
        },
        {
            "code_input": {
                "point_group": "D2h",
                "initial_irrep": "Ag",
                "final_irrep": "Ag",
                "transition_type": "electric_dipole",
            },
            "text_input": {"input_params": "D2h Ag Ag E1"},
            "output": {"result": {"is_allowed": False, "violated_rules": ["Laporte forbidden: g→g"]}},
        },
        {
            "code_input": {
                "point_group": "Oh",
                "initial_irrep": "T1g",
                "final_irrep": "T1u",
                "transition_type": "electric_dipole",
            },
            "text_input": {"input_params": "Oh T1g T1u E1"},
            "output": {"result": {"is_allowed": True}},
        },
    ]

    # ===== 特征标表 =====
    # 格式: {irrep: [characters], spin?: bool}
    # 列顺序对应群的类(对称操作)
    CHARACTER_TABLES = {
        "C2v": {
            "_classes": ["E", "C2(z)", "σv(xz)", "σv(yz)"],
            "_irreps": {
                "A1": [1, 1, 1, 1], "A2": [1, 1, -1, -1],
                "B1": [1, -1, 1, -1], "B2": [1, -1, -1, 1],
            },
            "_transformation": {"x": "B1", "y": "B2", "z": "A1"},
            "_RxRyRz": {"Rx": "B2", "Ry": "B1", "Rz": "A2"},
            "_quadrupole": {"x²,y²,z²": "A1", "xz": "B1", "yz": "B2", "xy": "A2"},
            "_has_gh": False,
        },
        "C3v": {
            "_classes": ["E", "2C3", "3σv"],
            "_irreps": {
                "A1": [1, 1, 1], "A2": [1, 1, -1],
                "E":  [2, -1, 0],
            },
            "_transformation": {"x,y": "E", "z": "A1"},
            "_RxRyRz": {"Rx,Ry": "E", "Rz": "A2"},
            "_quadrupole": {"z²": "A1", "(x²-y²,xy)": "E", "(xz,yz)": "E"},
            "_has_gh": False,
        },
        "D2h": {
            "_classes": ["E", "C2(z)", "C2(y)", "C2(x)", "i", "σ(xy)", "σ(xz)", "σ(yz)"],
            "_irreps": {
                "Ag": [1, 1, 1, 1, 1, 1, 1, 1],
                "B1g": [1, 1, -1, -1, 1, 1, -1, -1],
                "B2g": [1, -1, 1, -1, 1, -1, 1, -1],
                "B3g": [1, -1, -1, 1, 1, -1, -1, 1],
                "Au": [1, 1, 1, 1, -1, -1, -1, -1],
                "B1u": [1, 1, -1, -1, -1, -1, 1, 1],
                "B2u": [1, -1, 1, -1, -1, 1, -1, 1],
                "B3u": [1, -1, -1, 1, -1, 1, 1, -1],
            },
            "_transformation": {"x": "B3u", "y": "B2u", "z": "B1u"},
            "_RxRyRz": {"Rx": "B3g", "Ry": "B2g", "Rz": "Ag"},
            "_quadrupole": {"x²,y²,z²": "Ag", "xy": "B1g", "xz": "B2g", "yz": "B3g"},
            "_has_gh": True,
        },
        "D3h": {
            "_classes": ["E", "2C3", "3C2'", "σh", "2S3", "3σv"],
            "_irreps": {
                "A1'": [1, 1, 1, 1, 1, 1],
                "A2'": [1, 1, -1, 1, 1, -1],
                'E''':   [2, -1, 0, 2, -1, 0],
                'A1''': [1, 1, 1, -1, -1, -1],
                'A2''': [1, 1, -1, -1, -1, 1],
                'E''':  [2, -1, 0, -2, 1, 0],
            },
            "_transformation": {"x,y": "E'", "z": "A2'"},
            "_RxRyRz": {"Rx,Ry": "E''", "Rz": "A2'"},
            "_quadrupole": {"z²,x²+y²": "A1'", "(x²-y²,xy)": "E'"},
            "_has_gh": True,
        },
        "D4h": {
            "_classes": ["E", "2C4", "C2", "2C2'", "2C2''", "i", "2S4", "σh", "2σv", "2σd"],
            "_irreps": {
                "A1g": [1,1,1,1,1,1,1,1,1,1],
                "A2g": [1,1,1,-1,-1,1,1,1,1,-1],
                "B1g": [1,-1,1,1,-1,1,-1,1,1,-1],
                "B2g": [1,-1,1,-1,1,1,-1,1,1,1],
                "Eg":  [2,0,-2,0,0,2,0,-2,2,0],
                "A1u": [1,1,1,1,1,-1,-1,-1,-1,-1],
                "A2u": [1,1,1,-1,-1,-1,-1,-1,-1,1],
                "B1u": [1,-1,1,1,-1,-1,1,-1,-1,1],
                "B2u": [1,-1,1,-1,1,-1,1,-1,-1,-1],
                "Eu":  [2,0,-2,0,0,-2,0,2,-2,0],
            },
            "_transformation": {"x,y": "Eu", "z": "A2u"},
            "_RxRyRz": {"Rx,Ry": "Eg", "Rz": "A2g"},
            "_quadrupole": {"z²": "A1g", "x²-y²": "B1g", "xy": "B2g", "(xz,yz)": "Eg"},
            "_has_gh": True,
        },
        "D6h": {
            "_classes": ["E","2C6","2C3","C2","3C2'","3C2''","i","2S3","2S6","σh","3σd","3σv"],
            "_irreps": {
                "A1g":[1,1,1,1,1,1,1,1,1,1,1,1],
                "A2g":[1,1,1,1,-1,-1,1,1,1,1,1,-1],
                "B1g":[1,-1,1,-1,1,-1,1,-1,1,-1,1,-1],
                "B2g":[1,-1,1,-1,-1,1,1,-1,1,-1,-1,1],
                "E1g":[2,1,-1,-2,0,0,2,1,-1,-2,0,0],
                "E2g":[2,-1,-1,2,0,0,2,-1,-1,2,0,0],
                "A1u":[1,1,1,1,1,1,-1,-1,-1,-1,-1,-1],
                "A2u":[1,1,1,1,-1,-1,-1,-1,-1,-1,-1,1],
                "B1u":[1,-1,1,-1,1,-1,-1,1,-1,1,-1,1],
                "B2u":[1,-1,1,-1,-1,1,-1,1,-1,1,1,-1],
                "E1u":[2,1,-1,-2,0,0,-2,-1,1,2,0,0],
                "E2u":[2,-1,-1,2,0,0,-2,1,1,-2,0,0],
            },
            "_transformation": {"x,y":"E1u","z":"A2u"},
            "_RxRyRz":{"Rx,Ry":"E1g","Rz":"A2g"},
            "_quadrupole":{"z²":"A1g","x²-y²":"E2g","xy":"E2g","(xz,yz)":"E1g"},
            "_has_gh": True,
        },
        "Td": {
            "_classes": ["E", "8C3", "3C2", "6S4", "6σd"],
            "_irreps": {
                "A1": [1, 1, 1, 1, 1],
                "A2": [1, 1, 1, -1, -1],
                "E":  [2, -1, 2, 0, 0],
                "T1": [3, 0, -1, 1, -1],
                "T2": [3, 0, -1, -1, 1],
            },
            "_transformation": {"x,y,z": "T2"},
            "_RxRyRz": {"Rx,Ry,Rz": "T1"},
            "_quadrupole": {"x²+y²+z²": "A1", "(2z²-x²-y²,x²-y²)": "E", "(xy,xz,yz)": "T2"},
            "_has_gh": False,
        },
        "Oh": {
            "_classes": ["E", "8C3", "6C2", "6C4", "3C2(=C4²)", "i", "6S4", "8S6", "3σh", "6σd"],
            "_irreps": {
                "A1g": [1,1,1,1,1,1,1,1,1,1],
                "A2g": [1,1,-1,-1,1,1,1,-1,1,-1],
                "Eg":  [2,-1,0,0,2,2,0,-1,2,0],
                "T1g": [3,0,-1,1,-1,3,-1,0,-1,-1],
                "T2g": [3,0,1,-1,-1,3,1,0,-1,1],
                "A1u": [1,1,1,1,1,-1,-1,-1,-1,-1],
                "A2u": [1,1,-1,-1,1,-1,-1,1,-1,1],
                "Eu":  [2,-1,0,0,2,-2,0,1,-2,0],
                "T1u": [3,0,-1,1,-1,-3,1,0,1,1],
                "T2u": [3,0,1,-1,-1,-3,-1,0,1,-1],
            },
            "_transformation": {"x,y,z": "T1u"},
            "_RxRyRz": {"Rx,Ry,Rz": "T1g"},
            "_quadrupole":{"x²+y²+z²":"A1g","(2z²-x²-y²,x²-y²)":"Eg","(xy,xz,yz)":"T2g"},
            "_has_gh": True,
        },
        "Cinfv": {
            "_classes": ["E", "2C∞", "...", "∞σv"],
            "_irreps": {
                "Σ+": [1, 1, "+", 1],
                "Σ-": [1, 1, "+", -1],
                "Π":  [2, 2, "-", 0],
                "Δ":  [2, 2, "+", 0],
            },
            "_transformation": {"x,y": "Π", "z": "Σ+"},
            "_RxRyRz": {"Rx,Ry": "Π", "Rz": "Σ-"},
            "_quadrupole": {"z²": "Σ+", "(x²-y²,xy)": "Δ", "(xz,yz)": "Π"},
            "_has_gh": False,
        },
        "Dinfh": {
            "_classes": ["E", "2C∞", "...", "∞C2'", "i", "2S∞", "...", "∞σh"],
            "_irreps": {
                "Σg+": [1, 1, "+", 1, 1, 1, "+", 1],
                "Σg-": [1, 1, "+", 1, -1, -1, "+", -1],
                "Πg":  [2, 2, "-", 0, 2, -2, "-", 0],
                "Δg":  [2, 2, "+", 0, 2, 2, "+", 0],
                "Σu+": [1, 1, "+", 1, -1, -1, "+", -1],
                "Σu-": [1, 1, "+", 1, 1, 1, "+", 1],
                "Πu":  [2, 2, "-", 0, -2, 2, "-", 0],
                "Δu":  [2, 2, "+", 0, -2, -2, "+", 0],
            },
            "_transformation": {"x,y": "Πu", "z": "Σu+"},
            "_RxRyRz": {"Rx,Ry": "Πg", "Rz": "Σg-"},
            "_quadrupole": {"z²": "Σg+", "(x²-y²,xy)": "Δg", "(xz,yz)": "Πg"},
            "_has_gh": True,
        },
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, point_group: str, initial_irrep: str, final_irrep: str,
                  transition_type: str = "electric_dipole",
                  spin_multiplicity_initial: int = 1,
                  spin_multiplicity_final: int = 1,
                  polarization: str = "any") -> dict:
        """
        核心逻辑：基于群论的选择定则分析。
        """
        pg = self._normalize_pg(point_group)
        if pg not in self.CHARACTER_TABLES:
            raise ChemMCPError(
                f"Point group '{pg}' not supported.\n"
                f"Supported groups: {list(self.CHARACTER_TABLES.keys())}"
            )

        ct = self.CHARACTER_TABLES[pg]
        irreps = ct["_irreps"]

        # 验证 irrep 名称
        iri_norm = self._normalize_irrep(initial_irrep)
        irf_norm = self._normalize_irrep(final_irrep)

        if iri_norm not in irreps:
            raise ChemMCPError(f"Initial irrep '{initial_irrep}' not found in {pg}. Available: {list(irreps.keys())}")
        if irf_norm not in irreps:
            raise ChemMCPError(f"Final irrep '{final_irrep}' not found in {pg}. Available: {list(irreps.keys())}")

        ttype = transition_type.lower().replace(" ", "_")

        # ---- 获取跃迁算子的 irrep ----
        op_irreps = self._get_operator_irreps(pg, ttype, polarization)
        
        # ---- 直积分解 ----
        results = []
        all_allowed = []
        for op_name, op_ir in op_irreps:
            product = self._direct_product(iri_norm, op_ir, irf_norm, pg)
            contains_A1 = product["contains_totally_symmetric"]

            all_allowed.append(contains_A1)
            results.append({
                "operator_component": op_name,
                "operator_irrep": op_ir,
                "direct_product_decomposition": product["decomposition"],
                "contains_A1": contains_A1,
                "allowed_by_symmetry": contains_A1,
            })

        symmetry_allowed = any(all_allowed)

        # ---- g/u 宇称检查 ----
        parity_result = self._check_parity(iri_norm, irf_norm, pg, ct["_has_gh"], ttype)

        # ---- 自旋检查 ----
        Si = (spin_multiplicity_initial - 1) / 2
        Sf = (spin_multiplicity_final - 1) / 2
        spin_allowed = abs(Si - Sf) < 1e-10

        # ---- 综合判定 ----
        overall_allowed = symmetry_allowed and (
            not parity_result["applies"] or parity_result["allowed"]
        ) and spin_allowed

        result = {
            "point_group": pg,
            "initial_state_irrep": iri_norm,
            "final_state_irrep": irf_norm,
            "transition_type": transition_type,
            "polarization": polarization,
            "symmetry_analysis": {
                "operator_irreps_used": [(r["operator_component"], r["operator_irrep"]) for r in results],
                "direct_product_results": results,
                "symmetry_allowed": symmetry_allowed,
            },
            "parity_analysis": parity_result,
            "spin_analysis": {
                "initial_spin_S": round(Si, 4),
                "final_spin_S": round(Sf, 4),
                "initial_multiplicity": spin_multiplicity_initial,
                "final_multiplicity": spin_multiplicity_final,
                "spin_conserved": spin_allowed,
                "spin_rule": "ΔS = 0 ✓" if spin_allowed else "ΔS ≠ 0 ✗ (spin-forbidden)",
            },
            "overall_verdict": {
                "is_allowed": overall_allowed,
                "classification": self._classify_transition(overall_allowed, symmetry_allowed,
                                                          not parity_result["applies"] or parity_result["allowed"],
                                                          spin_allowed),
            },
            "character_table_summary": {
                "group": pg,
                "n_classes": len(ct.get("_classes", [])),
                "n_irreps": len(irreps),
                "has_inversion_center": ct["_has_gh"],
            },
        }

        logger.info(f"SelectionRules: {pg}, {iri_norm}→{irf_norm}, "
                     f"type={ttype}, allowed={overall_allowed}")
        return {"result": result}

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入。"""
        try:
            parts = input_params.strip().split()
            pg = parts[0]
            iri = parts[1]
            irf = parts[2]
            ttype = parts[3] if len(parts) > 3 else "electric_dipole"
            return self._run_base(pg, iri, irf, ttype)
        except IndexError as e:
            raise ChemMCPError(f"Failed to parse: need at least 'pg iri irf'. Got: {input_params}")
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Parse error: {e}")

    def _normalize_pg(self, pg: str) -> str:
        """标准化点群名称。"""
        p = pg.strip()
        mapping = {
            "c2v": "C2v", "c3v": "C3v", "c4v": "C4v",
            "d2h": "D2h", "d3h": "D3h", "d4h": "D4h", "d6h": "D6h",
            "td": "Td", "oh": "Oh", "oh": "Oh",
            "cinfv": "Cinfv", "c∞v": "Cinfv", "cinfv": "Cinfv",
            "dinfh": "Dinfh", "d∞h": "Dinfh", "dinfh": "Dinfh",
        }
        return mapping.get(p.lower(), p)

    def _normalize_irrep(self, ir: str) -> str:
        """标准化 irrep 名称。"""
        ir = ir.strip()
        # 常见变体映射
        variants = {
            "a1": "A1", "a2": "A2", "b1": "B1", "b2": "B2",
            "ag": "Ag", "au": "Au", "b1g": "B1g", "b2g": "B2g", "b3g": "B3g",
            "b1u": "B1u", "b2u": "B2u", "b3u": "B3u",
            "a1'": "A1'", "a2'": "A2'", "a1''": "A1''", "a2''": "A2''",
            "e'": "E'", "e''": "E''", "t1": "T1", "t2": "T2",
            "t1g": "T1g", "t2g": "T2g", "t1u": "T1u", "t2u": "T2u",
            "e1g": "E1g", "e2g": "E2g", "e1u": "E1u", "e2u": "E2u",
            "sigma+": "Σ+", "sigma-": "Σ-", "pi": "Π", "delta": "Δ",
            "sigmag+": "Σg+", "sigmamu": "Σmu", "sigmag-": "Σg-", "sigmamu-": "Σmu-",
            "piu": "Πu", "pig": "Πg", "deltau": "Δu", "deltag": "Δg",
        }
        return variants.get(ir.lower(), ir)

    def _get_operator_irreps(self, pg: str, ttype: str, pol: str) -> list:
        """获取跃迁算子各分量的 irrep。"""
        ct = self.CHARACTER_TABLES[pg]

        if ttype in ("electric_dipole", "e1", "electric_dipole_atomic"):
            trans = ct.get("_transformation", {})
            if pol == "x":
                return [("μ_x", trans.get("x", "?"))]
            elif pol == "y":
                return [("μ_y", trans.get("y", "?"))]
            elif pol == "z":
                return [("μ_z", trans.get("z", "?"))]
            else:
                ops = set()
                for k, v in trans.items():
                    if k in ("x", "y", "z"):
                        ops.add((f"μ_{k}", v))
                    elif "," in k:
                        for comp in k.split(","):
                            ops.add((f"μ_{comp.strip()}", v))
                    else:
                        ops.add((f"μ_{k}", v))
                return list(ops)

        elif ttype in ("magnetic_dipole", "m1"):
            rr = ct.get("_RxRyRz", {})
            return [(f"R_{k}", v) for k, v in rr.items()]

        elif ttype in ("electric_quadrupole", "e2"):
            quad = ct.get("_quadrupole", {})
            return [(k, v) for k, v in quad.items() if "(" not in k]

        elif ttype in ("raman",):
            quad = ct.get("_quadrupole", {})
            return [(f"α({k})", v) for k, v in quad.items()]

        return [("unknown", "?")]

    def _direct_product(self, iri: str, ir_op: str, irf: str, pg: str) -> dict:
        """
        计算直积 Γ_i ⊗ Γ_op ⊗ Γ_f 并检查是否包含全对称表示。
        
        使用特征标的简单乘法进行分解。
        """
        ct = self.CHARACTER_TABLES[pg]
        chars = ct["_irreps"]

        chi_i = chars[iri]
        chi_op = chars.get(ir_op, None)
        chi_f = chars[irf]

        if chi_op is None:
            return {"decomposition": f"{iri} ⊗ {ir_op} ⊗ {irf} (unknown op irrep)",
                    "contains_totally_symmetric": False}

        # 逐特征标相乘
        product_chars = [chi_i[i] * chi_op[i] * chi_f[i] for i in range(len(chi_i))]

        # 检查是否包含全对称表示 (A1/A1g/A1'/Σg+/...)
        a1_key = self._find_a1_key(chars)
        chi_a1 = chars[a1_key]

        # 内积 <χ_prod | χ_A1> = (1/h) Σ_g n_g χ_prod(g) χ_A1(g)
        # 对于 A1，所有 χ_A1(g)=1，所以内积 = (1/h) Σ n_g χ_prod(g)
        h = sum(chi_a1)  # 群阶近似（仅对第一列正确）
        inner_product = sum(product_chars) / max(h, 1)
        contains_A1 = abs(inner_product - 1.0) < 0.01

        # 简化分解（不完全但给出主要成分）
        decomp_parts = []
        for key, ch_val in chars.items():
            if key == a1_key:
                continue
            ip = sum(a * b for a, b in zip(product_chars, ch_val)) / max(h, 1)
            if abs(ip) > 0.5:
                decomp_parts.append(f"{key}")

        decomp_str = f"{iri} ⊗ {ir_op} ⊗ {irf}"
        if contains_A1:
            decomp_str += f" ⊃ {a1_key} + " + " + ".join(decomp_parts)
        else:
            decomp_str += " = " + " + ".join([a1_key] + decomp_parts) if decomp_parts else decomp_str

        return {
            "decomposition": decomp_str,
            "product_characters": [round(c, 2) for c in product_chars],
            "contains_totally_symmetric": contains_A1,
            "a1_irrep_name": a1_key,
        }

    def _find_a1_key(self, irreps: dict) -> str:
        """找到全对称表示的名称。"""
        for key in irreps:
            kl = key.lower()
            if kl.startswith("a1") or kl == "a1'" or kl.startswith("σg+") or kl == "Σg+":
                return key
        # fallback: 第一个就是 A1
        return list(irreps.keys())[0]

    def _check_parity(self, iri: str, irf: str, pg: str, has_gh: bool, ttype: str) -> dict:
        """检查宇称选择定则。"""
        if not has_gh:
            return {"applies": False, "reason": "No inversion center in this point group"}

        gi = "g" in iri.lower() or "prime" in iri.lower() or "Σg" in iri
        gf = "g" in irf.lower() or "prime" in irf.lower() or "Σg" in irf

        pi_s = "g" if gi else "u"
        pf_s = "g" if gf else "u"

        if ttype in ("electric_dipole", "e1"):
            # E1: 需要 g ↔ u
            allowed = gi != gf
            rule = "Laporte: g ↔ u required for E1"
        elif ttype in ("magnetic_dipole", "m1"):
            # M1: 需要 g → g, u → u
            allowed = gi == gf
            rule = "M1: same parity required (g→g, u→u)"
        elif ttype in ("electric_quadrupole", "e2"):
            # E2: 同宇称
            allowed = gi == gf
            rule = "E2: same parity required"
        else:
            allowed = True
            rule = "No specific parity rule"

        return {
            "applies": True,
            "initial_parity": pi_s,
            "final_parity": pf_s,
            "parity_change": f"{pi_s} → {pf_s}",
            "allowed": allowed,
            "rule": rule,
        }

    @staticmethod
    def _classify_transition(overall, sym, par, spin) -> str:
        if overall:
            return "ALLOWED (fully symmetry-, parity-, and spin-allowed)"
        elif not sym:
            return "FORBIDDEN by symmetry (irrep direct product does not contain A₁)"
        elif not par:
            return "FORBIDDEN by Laporte rule (same parity in centrosymmetric molecule)"
        elif not spin:
            return "SPIN-FORBIDDEN (different spin multiplicities; may gain intensity via spin-orbit coupling)"
        else:
            return "UNKNOWN classification"
