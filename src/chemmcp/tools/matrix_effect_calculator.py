import logging
import math
from typing import List, Optional

import numpy as np

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MatrixEffectCalculator(BaseTool):
    """
    基质效应定量评估工具。
    
    用于 LC-MS/MS 等分析方法验证中的基质效应评估（post-extraction addition法）。
    """
    __version__ = "0.1.0"
    name             = "MatrixEffectCalculator"
    func_name        = "calculate_matrix_effect"
    description      = "Quantitatively assess matrix effect (ME) for bioanalytical/LC-MS method validation using post-extraction spike method."
    implementation_description = "Calculates matrix effect as ME% = (A_post_extraction / A_neat_solvent - 1) × 100%. Also computes process efficiency (PE) and recovery (RE) when pre-extraction data provided."
    oss_dependencies = [
        ("NumPy", "https://numpy.org", "BSD"),
    ]
    services_and_software = []
    categories       = ["General"]
    tags             = ["Matrix Effect", "LC-MS", "Method Validation", "Bioanalysis"]
    required_envs    = []

    code_input_sig   = [
        ("neat_solvent_areas", "list", "N/A", "Peak areas of analyte in neat solvent (set A)."),
        ("post_extraction_areas", "list", "N/A", "Peak areas of analyte spiked into post-extracted blank matrix (set B)."),
        ("pre_extraction_areas_or_None", "list_or_None", "None", "Peak areas of analyte spiked before extraction (set C, optional for RE/PE)."),
        ("threshold_ion_suppression", "float", "-20.0", "Threshold (%) below which ion suppression is considered significant."),
        ("threshold_enhancement", "float", "20.0", "Threshold (%) above which ion enhancement is considered significant."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'neat_areas post_extr_areas [pre_extr_areas] [supp_thresh] [enh_thresh]'"),
    ]

    output_sig       = [
        ("result", "dict", "Dictionary with ME%, PE%, RE%, assessment, and detailed results per sample."),
    ]

    examples         = [
        {
            "code_input": {
                "neat_solvent_areas": [1000000, 980000, 1020000, 995000, 1010000],
                "post_extraction_areas": [890000, 870000, 910000, 885000, 900000],
                "pre_extraction_areas_or_None": [750000, 730000, 770000, 745000, 755000],
                "threshold_ion_suppression": -20.0,
                "threshold_enhancement": 20.0,
            },
            "text_input": {
                "input_params": "1000000,980000,1020000,995000,1010000 890000,870000,910000,885000,900000 750000,730000,770000,745000,755000",
            },
            "output": {
                "result": {
                    "mean_me_pct": "...",
                    "assessment": "...",
                    "mean_pe_pct": "...",
                    "mean_re_pct": "...",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        neat_solvent_areas: List[float],
        post_extraction_areas: List[float],
        pre_extraction_areas_or_None: Optional[List[float]] = None,
        threshold_ion_suppression: float = -20.0,
        threshold_enhancement: float = 20.0,
    ) -> dict:
        """核心逻辑：基质效应计算"""
        if len(neat_solvent_areas) != len(post_extraction_areas):
            raise ChemMCPError("Neat solvent and post-extraction area lists must have same length.")
        if len(neat_solvent_areas) < 3:
            raise ChemMCPError("At least 3 replicates are recommended.")

        a_neat = np.array(neat_solvent_areas, dtype=float)
        a_post = np.array(post_extraction_areas, dtype=float)
        n = len(a_neat)

        # ---- Matrix Effect (ME) ----
        # ME% = (A_post / A_neat - 1) × 100%
        me_values = ((a_post / a_neat) - 1.0) * 100.0
        mean_me = float(np.mean(me_values))
        sd_me = float(np.std(me_values, ddof=1))
        rsd_me = abs(mean_me) > 1e-10 and (sd_me / abs(mean_me) * 100) or 0.0

        # 评估
        if mean_me <= threshold_ion_suppression:
            assessment = f"Significant ion suppression (ME={mean_me:.1f}% ≤ {threshold_ion_suppression}%)"
        elif mean_me >= threshold_enhancement:
            assessment = f"Significant ion enhancement (ME={mean_me:.1f}% ≥ {threshold_enhancement}%)"
        else:
            assessment = f"Acceptable matrix effect (ME={mean_me:.1f}%)"

        result = {
            "n_replicates": n,
            "matrix_effect": {
                "mean_me_pct": round(mean_me, 4),
                "sd_me_pct": round(sd_me, 4),
                "rsd_me_pct": round(rsd_me, 4),
                "min_me_pct": round(float(np.min(me_values)), 4),
                "max_me_pct": round(float(np.max(me_values)), 4),
                "assessment": assessment,
            },
            "thresholds": {
                "suppression_threshold": threshold_ion_suppression,
                "enhancement_threshold": threshold_enhancement,
            },
            "sample_details": [],
        }

        # 每个样品详情
        for i in range(n):
            me_i = round(float(me_values[i]), 4)
            if me_i <= threshold_ion_suppression:
                effect = "suppression"
            elif me_i >= threshold_enhancement:
                effect = "enhancement"
            else:
                effect = "negligible"
            result["sample_details"].append({
                "sample": int(i + 1),
                "neat_area": round(float(a_neat[i]), 2),
                "post_extraction_area": round(float(a_post[i]), 2),
                "me_pct": me_i,
                "effect_type": effect,
            })

        # ---- Process Efficiency (PE) & Recovery (RE) ----
        if pre_extraction_areas_or_None is not None:
            a_pre = np.array(pre_extraction_areas_or_None, dtype=float)
            if len(a_pre) == n:
                # RE% = (A_pre / A_post) × 100%
                re_values = (a_pre / a_post) * 100.0
                # PE% = (A_pre / A_neat) × 100%
                pe_values = (a_pre / a_neat) * 100.0

                result["recovery"] = {
                    "mean_re_pct": round(float(np.mean(re_values)), 4),
                    "sd_re_pct": round(float(np.std(re_values, ddof=1)), 4),
                }
                result["process_efficiency"] = {
                    "mean_pe_pct": round(float(np.mean(pe_values)), 4),
                    "sd_pe_pct": round(float(np.std(pe_values, ddof=1)), 4),
                }
                # 更新详情
                for i in range(n):
                    result["sample_details"][i]["pre_extraction_area"] = round(float(a_pre[i]), 2)
                    result["sample_details"][i]["recovery_pct"] = round(float(re_values[i]), 4)
                    result["sample_details"][i]["process_efficiency_pct"] = round(float(pe_values[i]), 4)

        logger.info(f"Matrix effect evaluation: ME={mean_me:.2f}%, assessment={assessment}")
        return result

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            neat = [float(x) for x in parts[0].split(",")]
            post = [float(x) for x in parts[1].split(",")]
            pre = None
            supp_t = -20.0
            enh_t = 20.0
            if len(parts) > 2 and parts[2].lower() != "none":
                pre = [float(x) for x in parts[2].split(",")]
            if len(parts) > 3:
                supp_t = float(parts[3])
            if len(parts) > 4:
                enh_t = float(parts[4])
            return self._run_base(neat, post, pre, supp_t, enh_t)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'neat_areas post_areas [pre_areas] [supp_thresh] [enh_thresh]'")
