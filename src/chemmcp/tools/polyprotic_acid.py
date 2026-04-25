import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# 常见多元酸数据库
_POLYPROTIC_ACID_DB = {
    "carbonic acid":       {"formula": "H2CO3",   "pKa": [6.35, 10.33],     "n_protons": 2},
    "h2co3":               {"formula": "H2CO3",   "pKa": [6.35, 10.33],     "n_protons": 2},
    "phosphoric acid":     {"formula": "H3PO4",   "pKa": [2.15, 7.20, 12.35],"n_protons": 3},
    "h3po4":               {"formula": "H3PO4",   "pKa": [2.15, 7.20, 12.35],"n_protons": 3},
    "sulfuric acid":       {"formula": "H2SO4",   "pKa": [-3.0, 1.99],      "n_protons": 2},
    "h2so4":               {"formula": "H2SO4",   "pKa": [-3.0, 1.99],      "n_protons": 2},
    "oxalic acid":         {"formula": "H2C2O4",  "pKa": [1.25, 4.27],      "n_protons": 2},
    "h2c2o4":              {"formula": "H2C2O4",  "pKa": [1.25, 4.27],      "n_protons": 2},
    "citric acid":         {"formula": "C6H8O7",  "pKa": [3.13, 4.76, 6.40],"n_protons": 3},
    "hydrosulfuric acid":  {"formula": "H2S",     "pKa": [7.04, 19.0],      "n_protons": 2},
    "h2s":                 {"formula": "H2S",     "pKa": [7.04, 19.0],      "n_protons": 2},
    "chromic acid":        {"formula": "H2CrO4",  "pKa": [0.74, 6.49],      "n_protons": 2},
    "arsenic acid":        {"formula": "H3AsO4",  "pKa": [2.26, 6.76, 11.29],"n_protons": 3},
    "tartaric acid":       {"formula": "C4H6O6",  "pKa": [3.04, 4.37],      "n_protons": 2},
    "malic acid":          {"formula": "C4H6O5",  "pKa": [3.46, 5.11],      "n_protons": 2},
    "succinic acid":       {"formula": "C4H6O4",  "pKa": [4.21, 5.64],      "n_protons": 2},
    "phthalic acid":       {"formula": "C8H6O4",  "pKa": [2.89, 5.51],      "n_protons": 2},
    "boric acid (full)":   {"formula": "H3BO3",   "pKa": [9.24, 12.74, 13.80],"n_protons": 3},
}


@ChemMCPManager.register_tool
class PolyproticAcid(BaseTool):
    """
    多元酸分步解离计算。
    支持二元酸和三元酸，计算各物种浓度、pH 和分布分数（α 图数据）。
    """
    __version__ = "0.1.0"
    name = "PolyproticAcid"
    func_name = "polyprotic_acid_dissociation"
    description = "Calculate stepwise dissociation of polyprotic acids: pH, species concentrations, and fractional composition (α diagram data)."
    implementation_description = "Implements mass-balance and equilibrium equations for diprotic and triprotic acids to compute all species concentrations."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Polyprotic Acid", "Dissociation", "Species Distribution", "Equilibrium"]
    required_envs = []

    code_input_sig = [
        ("acid_name_or_pka_list", "str", "N/A", "Acid name from database OR comma-separated pKa values (e.g., '2.15,7.20,12.35')."),
        ("C0", "float", "0.1", "Initial total concentration of the polyprotic acid (mol/L)."),
        ("num_ph_points", "int", "50", "Number of pH points for fractional composition diagram."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'acid_name_or_pKas C0 [num_ph_points]'. Example: 'phosphoric acid 0.1 50'"),
    ]

    output_sig = [
        ("acid_name", "str", "Name of the acid used."),
        ("formula", "str", "Molecular formula."),
        ("pKa_values", "list", "List of pKa values."),
        ("n_protons", "int", "Number of dissociable protons."),
        ("ph", "float", "Calculated pH of the solution at C0."),
        ("species", "dict", "Concentrations of each species {species_name: conc_mol_L}."),
        ("alpha_diagram", "list", "Fractional composition data: list of {ph, species_fractions} for plotting."),
        ("dominant_species", "str", "The dominant species at the calculated pH."),
        ("explanation", "str", "Step-by-step explanation of the calculation."),
    ]

    examples = [
        {
            "code_input": {
                "acid_name_or_pka_list": "phosphoric acid",
                "C0": 0.1,
                "num_ph_points": 50,
            },
            "text_input": {
                "input_params": "phosphoric acid 0.1"
            },
            "output": {
                "acid_name": "phosphoric acid",
                "formula": "H3PO4",
                "pKa_values": [2.15, 7.20, 12.35],
                "n_protons": 3,
                "ph": 1.53,
                "species": {"H3PO4": 0.099, "H2PO4-": 0.001, "HPO4^2-": 0.0, "PO4^3-": 0.0},
                "alpha_diagram": [{"ph": 0.0, "fractions": {}}],
                "dominant_species": "H3PO4",
                "explanation": "Polyprotic acid dissociation.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.Kw = 1.0e-14

    def _run_base(self, acid_name_or_pka_list: str, C0: float = 0.1, num_ph_points: int = 50) -> dict:
        """核心逻辑：多元酸分步解离"""
        if C0 <= 0:
            raise ChemMCPError("Concentration must be positive.")

        input_str = str(acid_name_or_pka_list).strip()
        entry = self._resolve_input(input_str)
        pKa_list = entry["pKa"]
        n = entry["n_protons"]
        formula = entry["formula"]
        name = entry.get("name", input_str)

        Ka_list = [10 ** (-pKa) for pKa in pKa_list]

        # 计算溶液 pH
        ph = self._solve_ph(C0, Ka_list, n)

        # 计算各物种浓度
        h = 10 ** (-ph)
        species = self._calc_species(C0, h, Ka_list, n, name)

        # 分布分数图数据
        alpha_data = self._generate_alpha_diagram(Ka_list, n, num_ph_points, name)

        dominant = max(species.items(), key=lambda x: x[1])[0] if species else "unknown"

        explanation = (
            f"{name} ({formula}) is a {n}-protic acid with pKa values: "
            f"{', '.join(str(p) for p in pKa_list)}.\n"
            f"At C0 = {C0} M, pH = {ph:.2f}. Dominant species: {dominant}."
        )

        logger.info(f"Polyprotic acid {name}: pH={ph:.2f}, dominant={dominant}")

        return {
            "acid_name": name,
            "formula": formula,
            "pKa_values": pKa_list,
            "n_protons": n,
            "ph": round(ph, 4),
            "species": {k: round(v, 10) for k, v in species.items()},
            "alpha_diagram": alpha_data,
            "dominant_species": dominant,
            "explanation": explanation,
        }

    def _resolve_input(self, input_str: str) -> dict:
        key = input_str.lower()
        if key in _POLYPROTIC_ACID_DB:
            return {**_POLYPROTIC_ACID_DB[key], "name": input_str}
        matches = [k for k in _POLYPROTIC_ACID_DB if key in k or k in key]
        if len(matches) == 1:
            return {**_POLYPROTIC_ACID_DB[matches[0]], "name": matches[0]}
        elif len(matches) > 1:
            raise ChemMCPError(f"Multiple matches: {matches}. Be more specific.")
        try:
            pkas = [float(x.strip()) for x in input_str.split(",")]
            if len(pkas) < 2:
                raise ValueError("Need at least 2 pKa values for polyprotic acid.")
            return {"name": f"custom ({len(pkas)}-protic)", "formula": f"H{len(pkas)+1}A",
                    "pKa": pkas, "n_protons": len(pkas)}
        except ValueError:
            available = sorted(_POLYPROTIC_ACID_DB.keys())
            raise ChemMCPError(f"Unknown acid '{input_str}'. Available: {available}\nOr provide comma-separated pKa values.")

    def _compute_D_and_alphas(self, h_val: float, Ka_list: list, n: int) -> tuple:
        """Compute denominator D and all alpha fractions."""
        D = h_val ** n
        terms = [D]  # terms[i] = numerator for α_i
        for i in range(1, n + 1):
            term = 1.0
            for j in range(i):
                term *= Ka_list[j]
            term *= h_val ** (n - i)
            terms.append(term)
            D += term
        D = max(D, 1e-300)
        alphas = [t / D for t in terms]
        return alphas

    def _solve_ph(self, C0: float, Ka_list: list, n: int) -> float:
        """求解 pH：先粗扫描找符号变化区间，再二分法精确定位"""
        def f(h_val):
            alphas = self._compute_D_and_alphas(h_val, Ka_list, n)
            oh = self.Kw / h_val
            anion_charge = sum(i * alphas[i] * C0 for i in range(1, n + 1))
            return h_val - oh - anion_charge

        # 粗扫描：在 pH 0-14 范围内找 f(h) 的符号变化
        ph_low, ph_high = 0.0, 0.0
        prev_val = None
        for ph_test in [i * 0.1 for i in range(141)]:
            h_test = 10 ** (-ph_test)
            val = f(h_test)
            if prev_val is not None and prev_val * val < 0:
                ph_low, ph_high = ph_test - 0.1, ph_test
                break
            prev_val = val

        if ph_high == 0.0:
            # 没有找到符号变化，回退到近似公式（用第一级解离）
            Ka1 = Ka_list[0]
            disc = Ka1**2 + 4*Ka1*C0
            h = (-Ka1 + math.sqrt(disc)) / 2
            return -math.log10(max(h, 1e-14))

        # 在找到的区间内用二分法精炼
        for _ in range(100):
            mid = (ph_low + ph_high) / 2
            h_mid = 10 ** (-mid)
            val = f(h_mid)
            if val > 0:
                ph_low = mid
            else:
                ph_high = mid
            if ph_high - ph_low < 1e-10:
                break
        return (ph_low + ph_high) / 2

    def _calc_species(self, C0: float, h: float, Ka_list: list, n: int, name: str) -> dict:
        """计算各物种浓度"""
        alphas = self._compute_D_and_alphas(h, Ka_list, n)
        species_names = self._get_species_names(name, n)
        species = {}
        for i in range(n + 1):
            sname = species_names[i] if i < len(species_names) else f"A^{i-n:-d}"
            species[sname] = alphas[i] * C0
        return species

    def _get_species_names(self, name: str, n: int) -> list:
        base = name.split("(")[0].strip().title().replace(" ", "")
        if n == 2:
            if "carbonic" in name.lower() or "h2co3" in name.lower():
                return ["H2CO3", "HCO3-", "CO3^2-"]
            elif "oxalic" in name.lower() or "h2c2o4" in name.lower():
                return ["H2C2O4", "HC2O4-", "C2O4^2-"]
            elif "sulfuric" in name.lower() or "h2so4" in name.lower():
                return ["H2SO4", "HSO4-", "SO4^2-"]
            return [f"H2{base}", f"H{base}-", f"{base}2-"]
        elif n == 3:
            if "phosphoric" in name.lower() or "h3po4" in name.lower():
                return ["H3PO4", "H2PO4-", "HPO4^2-", "PO4^3-"]
            elif "citric" in name.lower():
                return ["H3Cit", "H2Cit-", "HCit^2-", "Cit^3-"]
            elif "arsenic" in name.lower():
                return ["H3AsO4", "H2AsO4-", "HAsO4^2-", "AsO4^3-"]
            return [f"H3{base}", f"H2{base}-", f"H{base}2-", f"{base}3-"]
        return [f"H{n-i+1}A{'+' if n-i>0 else ''}" for i in range(n+1)]

    def _generate_alpha_diagram(self, Ka_list: list, n: int, num_points: int, name: str) -> list:
        species_names = self._get_species_names(name, n)
        pka_max = max([-math.log10(k) for k in Ka_list])
        ph_min, ph_max = 0.0, min(14.0, pka_max + 2.0)
        data = []
        for i in range(num_points):
            ph = ph_min + (ph_max - ph_min) * i / (num_points - 1) if num_points > 1 else 7.0
            h = 10 ** (-ph)
            alphas = self._compute_D_and_alphas(h, Ka_list, n)
            fractions = {}
            for idx in range(len(species_names)):
                sname = species_names[idx]
                fractions[sname] = round(alphas[idx], 10)
            data.append({"ph": round(ph, 2), "fractions": fractions})
        return data

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            if len(parts) < 1:
                raise ValueError("Need acid name or pKa values.")
            acid_input = parts[0]
            C0 = float(parts[1]) if len(parts) > 1 else 0.1
            npts = int(parts[2]) if len(parts) > 2 else 50
            return self._run_base(acid_input, C0, npts)
        except (ValueError, IndexError) as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'acid_name/pKas [C0] [npts]'")
