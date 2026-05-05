import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Derivatization reagent database
DERIVATIZATION_REAGENTS = [
    {
        "reagent": "BSTFA + TMCS (1%)",
        "full_name": "N,O-Bis(trimethylsilyl)trifluoroacetamide + Trimethylchlorosilane",
        "target_functional_groups": ["-OH", "-COOH", "-NH2", "-SH", "-CONH2"],
        "target_analytes": ["acids", "alcohols", "phenols", "amines", "sugars", "steroids", "fatty_acids", "metabolites"],
        "detection_methods": ["GC-MS", "GC-FID", "GC-ECD"],
        "reaction_type": "Silylation (TMS derivative)",
        "reaction_temp_c": "60-80",
        "reaction_time_min": "15-30",
        "solvent": "Pyridine, acetonitrile, DMF, or ethyl acetate",
        "mechanism": "Replaces active H with -Si(CH3)3 group; increases volatility and thermal stability.",
        "pros": ["Universal for most polar groups", "Fast reaction", "Volatile byproducts", "Well-established"],
        "cons": ["Moisture sensitive", "Pyridine odor", "Some derivatives may be unstable"],
        "notes": "Most popular silylation reagent. TMCS acts as catalyst.",
    },
    {
        "reagent": "MSTFA",
        "full_name": "N-Methyl-N-(trimethylsilyl)trifluoroacetamide",
        "target_functional_groups": ["-OH", "-COOH", "-NH2", "-SH"],
        "target_analytes": ["acids", "alcohols", "amines", "steroids", "drugs_of_abuse", "metabolites"],
        "detection_methods": ["GC-MS", "GC-MS/MS"],
        "reaction_type": "Silylation (TMS derivative)",
        "reaction_temp_c": "60-80",
        "reaction_time_min": "10-20",
        "solvent": "Pyridine-free options available; acetonitrile, ethyl acetate",
        "mechanism": "Similar to BSTFA but with methyl group instead of O; fewer side products.",
        "pros": ["Pyridine-free formulations available", "Fewer byproducts than BSTFA", "Excellent for biological samples"],
        "cons": ["More expensive than BSTFA", "Still moisture sensitive"],
        "notes": "Preferred for GC-MS/MS of drugs/metabolites in biological matrices.",
    },
    {
        "reagent": "PFBBr (Pentafluorobenzyl bromide)",
        "full_name": "Pentafluorobromobenzene (α-Bromo-2,3,4,5,6-pentafluorotoluene)",
        "target_functional_groups": ["-COOH", "-OH (phenolic)", "-SH"],
        "target_analytes": ["carboxylic_acids", "phenols", "thiols", "halogenated_acids"],
        "detection_methods": ["GC-ECD", "GC-NCI-MS", "GC-ECNI-MS"],
        "reaction_type": "Alkylation (PFB ester/ether)",
        "reaction_temp_c": "50-70",
        "reaction_time_min": "30-60",
        "solvent": "Acetone, acetonitrile, DMF, with K2CO3 as base catalyst",
        "mechanism": "Introduces pentafluorobenzyl group — highly electron-capturing for ECD/NCI detection.",
        "pros": ["Excellent sensitivity in ECD/NCI mode", "Stable derivatives", "Low detection limits possible"],
        "cons": ["Specific to acidic/phenolic groups", "Longer reaction time", "Requires basic conditions"],
        "notes": "Gold standard for carboxylic acid analysis by GC-ECD or GC-NCI-MS.",
    },
    {
        "reagent": "DNPH (2,4-Dinitrophenylhydrazine)",
        "full_name": "2,4-Dinitrophenylhydrazine",
        "target_functional_groups": ["C=O (aldehyde)", "C=O (ketone)"],
        "target_analytes": ["formaldehyde", "acetaldehyde", "other_aldehydes", "ketones", "carbonyl_compounds"],
        "detection_methods": ["HPLC-UV", "LC-MS", "UV-Vis spectrophotometry"],
        "reaction_type": "Hydrazone formation",
        "reaction_temp_c": "25-40 (room temp)",
        "reaction_time_min": "30-120",
        "solvent": "Acidified acetonitrile or methanol (with acid catalyst)",
        "mechanism": "Forms colored hydrazone adducts with carbonyl compounds; strong UV absorption (~360 nm).",
        "pros": ["Simple reaction", "Strong UV chromophore", "EPA Method compliant for aldehydes/ketones"],
        "cons": ["Only for carbonyls", "Can have interference from other DNPH-reactive species", "Slow for some ketones"],
        "notes": "Standard method for formaldehyde and other carbonyl analysis in air/water samples.",
    },
    {
        "reagent": "Dansyl Chloride (Dns-Cl)",
        "full_name": "5-(Dimethylamino)naphthalene-1-sulfonyl chloride",
        "target_functional_groups": ["-NH2", "-NH-", "-OH (phenolic)", "-SH"],
        "target_analytes": ["primary_secondary_amines", "phenols", "thiols", "amino_acids", "biogenic_amines"],
        "detection_methods": ["HPLC-FLD", "LC-MS/MS", "HPLC-UV"],
        "reaction_type": "Sulfonylation (dansyl derivative)",
        "reaction_temp_c": "40-60",
        "reaction_time_min": "20-45",
        "solvent": "Acetone, acetonitrile, buffer pH 9-10",
        "mechanism": "Introduces fluorescent dansyl group enabling fluorescence detection at very low levels.",
        "pros": ["Excellent fluorescence properties", "Very high sensitivity (fmol level)", "Works for amines AND phenols"],
        "cons": ["Light sensitive", "Requires alkaline conditions", "Derivative stability varies"],
        "notes": "Go-to derivatization for biogenic amine and trace amine analysis by HPLC-FLD.",
    },
    {
        "reagent": "FMOC-Cl (Fluorenylmethyloxycarbonyl chloride)",
        "full_name": "9-Fluorenylmethyl chloroformate",
        "target_functional_groups": ["-NH2", "-NH- (secondary)"],
        "target_analytes": ["amino_acids", "amines", "peptides", "isocyanates"],
        "detection_methods": ["HPLC-FLD", "LC-MS/MS", "HPLC-UV"],
        "reaction_type": "Carbamylation (FMOC derivative)",
        "reaction_temp_c": "Room temperature",
        "reaction_time_min": "5-15 (very fast)",
        "solvent": "Acetonitrile/buffer (pH > 8), borate buffer common",
        "mechanism": "Introduces fluorenyl group — strongly fluorescent and UV-active.",
        "pros": ["Very fast reaction (minutes)", "Highly sensitive", "Excellent for amino acids", "Room temperature operation"],
        "cons": ["Excess reagent must be removed (interferes with chromatography)", "Hydrolyzes in water", "Can give multiple derivatives"],
        "notes": "Most popular for amino acid analysis. Must remove excess FMOC (extraction or secondary reaction).",
    },
    {
        "reagent": "OPA + Thiol (o-Phthalaldehyde)",
        "full_name": "o-Phthalaldehyde with thiol (e.g., NAC, MCE, mercaptoethanol)",
        "target_functional_groups": ["-NH2 (primary amine only)"],
        "target_analytes": ["primary_amines", "amino_acids", "biogenic_amines", "peptides"],
        "detection_methods": ["HPLC-FLD", "CE-LIF"],
        "reaction_type": "Condensation (isoindole formation)",
        "reaction_temp_c": "Room temperature",
        "reaction_time_min": "1-3 (very fast, pre-column) / instant (post-column)",
        "solvent": "Borate buffer pH 10.4 + thiol co-reagent",
        "mechanism": "OPA + primary amine + thiol → highly fluorescent isoindole product.",
        "pros": ["Extremely fast", "No secondary amines react (selective)", "Post-column compatible", "Widely used"],
        "cons": ["Primary amines only", "Derivatives are unstable (must analyze quickly)", "Thiol odor"],
        "notes": "Classic pre-column/post-column derivatization for amino acids. Derivatives degrade within minutes — use auto-sampler.",
    },
    {
        "reagent": "Marfey's Reagent (FDAA)",
        "full_name": "1-Fluoro-2,4-dinitrophenyl-5-L-alanine amide (FDAA)",
        "target_functional_groups": ["-NH2 (of amino acids)"],
        "target_analytes": ["amino_acids", "amines_for_chirality_analysis"],
        "detection_methods": ["HPLC-UV", "LC-MS"],
        "reaction_type": "Amide formation (diastereomer)",
        "reaction_temp_c": "40",
        "reaction_time_min": "60-120",
        "solvent": "1M NaHCO3 / acetone mixture",
        "mechanism": "Forms diastereomeric amides separable on reversed-phase HPLC for chiral resolution.",
        "pros": ["Gold standard for amino acid chirality determination", "UV active", "Relatively simple"],
        "cons": ["Long reaction time", "Limited to amino acids/amines", "Expensive reagent"],
        "notes": "Essential tool for determining D/L configuration of amino acids without chiral column.",
    },
    {
        "reagent": "TMO (Trimethyloxonium tetrafluoroborate)",
        "full_name": "Trimethyloxonium tetrafluoroborate (Meerwein's reagent)",
        "target_functional_groups": ["-COOH", "-PO4", "-SO4", "-OH (acidic)"],
        "target_analytes": ["carboxylic_acids", "phosphates", "sulfonates", "organic_acids"],
        "detection_methods": ["GC-MS", "GC-FID"],
        "reaction_type": "Methylation (methyl ester/ether)",
        "reaction_temp_c": "Room temperature to 50",
        "reaction_time_min": "15-30",
        "solvent": "DCM, dichloromethane, with base if needed",
        "mechanism": "Powerful methylating agent — converts acids to methyl esters under mild conditions.",
        "pros": ["Mild conditions (room temp)", "Fast", "No aqueous workup needed", "Quantitative"],
        "cons": ["Moisture sensitive", "Toxic/handling precautions needed", "Limited commercial availability"],
        "notes": "Alternative to diazomethane (much safer). Excellent for organic acid profiling.",
    },
    {
        "reagent": "BF3-Methanol",
        "full_name": "Boron trifluoride in methanol (14% w/v)",
        "target_functional_groups": ["-COOH (especially fatty acids)"],
        "target_analytes": ["fatty_acids", "fatty_acid_methyl_esters_FAMEs", "lipids"],
        "detection_methods": ["GC-FID", "GC-MS"],
        "reaction_type": "Transesterification / Esterification (methyl ester)",
        "reaction_temp_c": "70-100 (reflux)",
        "reaction_time_min": "30-90",
        "solvent": "BF3-methanol reagent solution (14%), sometimes with toluene cosolvent",
        "mechanism": "BF3 catalyzes conversion of fatty acids (free or glyceride-bound) to FAMEs.",
        "pros": ["Official AOAC/AOCS method", "Handles free fatty acids AND triglycerides", "Robust and reproducible"],
        "cons": ["High temperature required", "Corrosive reagent", "Not suitable for thermolabile analytes"],
        "notes": "Standard method for fatty acid analysis (FAME preparation) in food/oil chemistry.",
    },
]


@ChemMCPManager.register_tool
class DerivatizationReagentSelector(BaseTool):
    """
    衍生化试剂选择：根据目标分析物和检测方法推荐合适的衍生化试剂。
    """
    __version__ = "0.1.0"
    name = "DerivatizationReagentSelector"
    func_name = "select_derivatization_reagent"
    description = "Select the optimal derivatization reagent based on target functional groups, analyte class, and detection method."
    implementation_description = "Matches target analyte functional groups against a comprehensive derivatization reagent database with reaction conditions and protocols."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Derivatization", "Sample Preparation", "GC-MS", "HPLC", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("functional_group", "str", "N/A", "Target functional group to derivatize (e.g., '-COOH', '-OH', '-NH2', 'C=O', 'mixed')."),
        ("analyte_class", "str", "general", "Class of analytes: 'acids', 'alcohols', 'amines', 'amino_acids', 'carbonyls', 'fatty_acids', 'pesticides', 'drugs', 'general'."),
        ("detection_method", "str", "GC-MS", "Detection method after derivatization: 'GC-MS', 'GC-FID', 'GC-ECD', 'HPLC-UV', 'HPLC-FLD', 'LC-MS', 'LC-MS/MS'."),
        ("priority", "str", "balanced", "Priority: 'sensitivity', 'speed', 'simplicity', 'stability', or 'balanced'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'functional_group [analyte_class] [detection_method] [priority]'"),
    ]

    output_sig = [
        ("recommended_reagent", "dict", "Best derivatization reagent with full details."),
        ("protocol", "list", "Step-by-step derivization protocol."),
        ("rationale", "str", "Why this reagent was recommended."),
        ("alternatives", "list", "Other suitable reagents ranked by relevance."),
    ]

    examples = [
        {
            "code_input": {
                "functional_group": "-COOH",
                "analyte_class": "fatty_acids",
                "detection_method": "GC-FID",
            },
            "text_input": {
                "input_params": "-COOH fatty_acids GC-FID",
            },
            "output": {
                "recommended_reagent": {"reagent": "BF3-Methanol"},
                "rationale": "BF3-methanol is the official AOAC method for converting fatty acids to FAMEs for GC-FID analysis.",
            },
        },
        {
            "code_input": {
                "functional_group": "-NH2",
                "analyte_class": "amino_acids",
                "detection_method": "HPLC-FLD",
            },
            "text_input": {
                "input_params": "-NH2 amino_acids HPLC-FLD",
            },
            "output": {
                "recommended_reagent": {"reagent": "FMOC-Cl"},
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        functional_group: str,
        analyte_class: str = "general",
        detection_method: str = "GC-MS",
        priority: str = "balanced",
    ) -> dict:
        """Core logic: select best derivatization reagent."""
        fg = functional_group.upper().strip()
        ac = analyte_class.lower().strip()
        det = detection_method.upper().strip()
        pri = priority.lower().strip()

        # Score each reagent
        scored = []
        for r in DERIVATIZATION_REAGENTS:
            score = self._score_reagent(r, fg, ac, det, pri)
            scored.append((score, r))

        scored.sort(key=lambda x: x[0], reverse=True)
        best_score, best = scored[0]

        if best_score < 0.3:
            raise ChemMCPError(
                f"No suitable derivatization reagent found for functional group '{fg}' with {det} detection. "
                f"Available functional groups covered: -COOH, -OH, -NH2, -SH, C=O, -PO4, -SO4"
            )

        protocol = self._build_protocol(best)
        rationale = self._build_rationale(best, fg, ac, det)

        alternatives = []
        for score, r in scored[1:5]:
            if score >= 0.2:
                alternatives.append({
                    "reagent": r["reagent"],
                    "full_name": r["full_name"],
                    "score": round(score, 2),
                    "best_for": r.get("notes", ""),
                })

        logger.info(f"Derivatization reagent selected: {best['reagent']} for {fg} / {det}")
        return {
            "recommended_reagent": {k: v for k, v in best.items()},
            "protocol": protocol,
            "rationale": rationale,
            "alternatives": alternatives,
        }

    def _score_reagent(self, r, fg, ac, det, pri):
        score = 0.0

        # Functional group match (+4 max)
        fgs_upper = [g.upper() for g in r["target_functional_groups"]]
        if fg == "MIXED":
            if len(r["target_functional_groups"]) >= 3:
                score += 3.0
            else:
                score += 1.0
        elif fg in fgs_upper:
            score += 4.0
        elif any(fg.replace("-", "") in g.replace("-", "") for g in fgs_upper):
            score += 2.0

        # Analyte class match (+2 max)
        acs_lower = [a.lower() for a in r["target_analytes"]]
        if ac in acs_lower:
            score += 2.0
        elif any(ac.replace("_", "") in a.replace("_", "") for a in acs_lower):
            score += 1.0

        # Detection method match (+2 max)
        dets_upper = [d.upper() for d in r["detection_methods"]]
        if det in dets_upper:
            score += 2.0
        elif any(det.replace("-", "") in d.replace("-", "") for d in dets_upper):
            score += 1.0

        # Priority weighting
        if pri == "sensitivity" and any(w in r.get("notes", "").lower() for w in ["sensitive", "low detection"]):
            score += 1.0
        elif pri == "speed" and int(r.get("reaction_time_min", "30").split("-")[0]) <= 15:
            score += 1.0
        elif pri == "simplicity" and "simple" in r.get("notes", "").lower() or "room temp" in r.get("reaction_temp_c", "").lower():
            score += 1.0
        elif pri == "stability" and "stable" in r.get("notes", "").lower():
            score += 1.0

        return score

    def _build_protocol(self, r):
        return [
            f"1. Prepare sample: Ensure sample is dry (remove water if using silylation reagents).",
            f"2. Reagent preparation: Dissolve/dilute {r['reagent']} in {r['solvent']}. Typical concentration: 1-10 mg/mL.",
            f"3. Reaction: Mix sample with reagent (typically 2-10x molar excess). Incubate at {r['reaction_temp_c']}°C for {r['reaction_time_min']} min.",
            f"4. Quenching (if needed): Remove excess reagent or quench reaction per specific protocol.",
            f"5. Post-derivatization: Dilute with mobile phase or solvent compatible with injection.",
            f"6. Analysis: Inject into {', '.join(r['detection_methods'])}. Store derivatives appropriately (check stability).",
        ]

    def _build_rationale(self, r, fg, ac, det):
        parts = [f"{r['reagent']} selected:"]
        if fg.upper() in [g.upper() for g in r["target_functional_groups"]]:
            parts.append(f"It targets {fg} functional groups via {r['reaction_type']}.")
        if det.upper() in [d.upper() for d in r["detection_methods"]]:
            parts.append(f"It is optimized for {det} detection.")
        if ac.lower() in [a.lower() for a in r["target_analytes"]]:
            parts.append(f"Well-established for {ac} analysis.")
        return " ".join(parts)

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            fg = parts[0]
            ac = parts[1] if len(parts) > 1 else "general"
            det = parts[2] if len(parts) > 2 else "GC-MS"
            pri = parts[3] if len(parts) > 3 else "balanced"
            return self._run_base(fg, ac, det, pri)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'functional_group [analyte_class] [detection] [priority]'")
