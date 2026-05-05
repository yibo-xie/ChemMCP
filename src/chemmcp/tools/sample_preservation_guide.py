import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Sample preservation database
PRESERVATION_DATA = {
    "water_inorganic": {
        "name": "Water (Inorganic Analysis)",
        "recommended_container": "HDPE bottle (pre-cleaned)",
        "temperature": "4°C (refrigerated)",
        "max_holding_time": {
            "default": "28 days",
            "metals": "6 months with HNO3 pH<2",
            "nutrients": "48 h (refrigerated), 28 days (frozen)",
            "anions": "28-48 h",
            "cyanide": "14 days",
            "Cr(VI)": "24 h",
            "mercury": "28 days (oxidizing conditions)",
        },
        "preservatives": [
            {"analyte": "Metals (total)", "additive": "HNO3 to pH < 2", "reason": "Prevents adsorption to container walls and precipitation"},
            {"analyte": "Nutrients (N, P)", "additive": "H2SO4 to pH < 2, or freeze immediately", "reason": "Inhibits microbial activity"},
            {"analyte": "Cyanide", "additive": "NaOH to pH > 12, refrigerate", "reason": "Prevents loss as HCN gas"},
            {"analyte": "Sulfide", "additive": "Zn acetate + NaOH (pH > 9)", "reason": "Precipitates as ZnS, prevents volatilization"},
            {"analyte": "Phenols", "additive": "H2SO4 to pH < 2, refrigerate 4°C", "reason": "Prevents biodegradation"},
            {"analyte": "Oil & Grease", "additive": "H2SO4 to pH < 2, refrigerate 4°C", "reason": "Prevents biodegradation"},
            {"analyte": "Chromium VI", "additive": "Cool to 4°C, do NOT acidify!", "reason": "Acid converts Cr(VI) to Cr(III)"},
            {"analyte": "DO (Dissolved Oxygen)", "additive": "Fix immediately in field (Winkler reagents)", "reason": "Biological activity continues after sampling"},
        ],
        "notes": "Use containers rinsed with sample (triplicate rinse) before filling. Fill completely to minimize headspace.",
    },
    "water_organic": {
        "name": "Water (Organic / VOC Analysis)",
        "recommended_container": "Glass vial with PTFE-lined septum (zero headspace for VOCs)",
        "temperature": "4°C (VOCs: cool, no freezing)",
        "max_holding_time": {
            "voc": "14 days (EPA Method 524)", 
            "svoc": "40 days refrigerated",
            "pesticides": "7-14 days (check method-specific)",
            "pcb": "40 days refrigerated (7 days if extracted)",
            "pah": "14 days refrigerated, 30 days frozen",
        },
        "preservatives": [
            {"analyte": "VOCs", "additive": "HCl to pH < 2 (for some methods), refrigerate 4°C, zero headspace", "reason": "Inhibits biodegradation; prevent volatilization"},
            {"analyte": "SVOCs/pesticides", "additive": "Refrigerate at 4°C, protect from light", "reason": "Slow degradation; photolysis prevention"},
            {"analyte": "Herbicides (acidic)", "additive": "H2SO4 to pH ~2, refrigerate", "reason": "Stabilizes acidic herbicides"},
        ],
        "notes": "For VOCs: fill container completely with no bubbles. Do NOT filter before preservation. Use field blanks.",
    },
    "soil_sediment": {
        "name": "Soil / Sediment",
        "recommended_container": "Wide-mouth glass jar (amber) or pre-cleaned HDPE container",
        "temperature": "4°C (refrigerated); -20°C for long-term/volatile analytes",
        "max_holding_time": {
            "default": "6 months (frozen)",
            "voc": "14 days (4°C), 6 months (-20°C)",
            "svoc": "14 days (4°C), 6 months (-20°C)",
            "pesticides": "14-30 days (4°C), 1+ year (-20°C)",
            "metals": "Indefinite (air-dried and stored properly)",
        },
        "preservatives": [
            {"analyte": "VOCs", "additive": "Field-extract immediately or freeze at -20°C within hours", "reason": "Volatile loss from soil matrix"},
            {"analyte": "Mercury", "additive": "Store at 4°C, analyze within 28 days; add preservative vial if needed", "reason": "Prevents Hg loss/transformations"},
            {"analyte": "Sulfides/sulfates", "additive": "Refrigerate, minimize exposure to air", "reason": "Oxidation potential"},
        ],
        "notes": "Composite samples should be homogenized before subsampling. Document exact collection location and depth.",
    },
    "biological_fluid": {
        "name": "Biological Fluids (Blood, Plasma, Serum, Urine)",
        "recommended_container": "Polypropylene tube (blood: with anticoagulant if plasma)",
        "temperature": "-20°C (short-term); -80°C (long-term >1 month)",
        "max_holding_time": {
            "plasma_serum": "6 months at -80°C",
            "urine": "30 days at -20°C",
            "whole_blood": "Separate plasma/serum within 2 h of collection",
            "unstable_analytes": "Days to weeks depending on compound stability",
        },
        "preservatives": [
            {"analyte": "Plasma/Serum general", "additive": "Centrifuge within 30 min, aliquot, freeze at -80°C", "reason": "Prevents enzymatic degradation and protein changes"},
            {"analyte": "Urine general", "additive": "Refrigerate during collection, aliquot, freeze at -20°C", "reason": "Bacterial growth prevention"},
            {"analyte": "Drugs of abuse (urine)", "additive": "Adjust to pH 5-9, add sodium azide (0.1% w/v) if storing >48h", "reason": "Prevents bacterial degradation of drugs"},
            {"analyte": "Vitamins (light-sensitive)", "additive": "Protect from light (amber tubes), add antioxidant (ascorbic acid/BHT)", "reason": "Photodegradation and oxidation prevention"},
            {"analyte": "Glucose (blood)", "additive": "Fluoride/oxalate tube", "reason": "Inhibits glycolysis for up to 3 days"},
        ],
        "notes": "Avoid repeated freeze-thaw cycles (aliquot into single-use volumes). For enzyme analysis, keep on ice and process within 1 hour.",
    },
    "food_agricultural": {
        "name": "Food / Agricultural Products",
        "recommended_container": "Food-grade polyethylene bag or glass jar; frozen: freezer-safe container",
        "temperature": "Fresh: 4°C; Frozen: -20°C; Long-term: -80°C",
        "max_holding_time": {
            "fresh_produce": "Analyze within 24-48 h (store at 4°C)",
            "frozen_food": "6-12 months at -20°C",
            "pesticide_residue": "Analyze ASAP; freeze if delay >24h; stable months at -20°C",
            "mycotoxins": "Stable weeks-months at -20°C; years at -80°C",
            "nutrients_vitamins": "Highly variable — analyze ASAP, protect from light/oxygen",
        },
        "preservatives": [
            {"analyte": "Pesticides (non-fatty foods)", "additive": "Freeze at -20°C within 24 h of collection", "reason": "Enzymatic degradation slows dramatically when frozen"},
            {"analyte": "Pesticides (fatty foods)", "additive": "Freeze at -20°C, protect from light oxidation", "reason": "Lipid oxidation can affect pesticide stability"},
            {"analyte": "Vitamin C", "additive": "Add metaphosphoric acid (stabilizer), store in amber container, freeze", "reason": "Rapidly oxidizes; acid and cold slow degradation"},
            {"analyte": "Fatty acids (for profile)", "additive": "Flush with N2/Ar gas, add BHT (antioxidant), store at -80°C", "reason": "Prevent lipid peroxidation"},
        ],
        "notes": "Homogenize representative subsample before analysis. Document storage conditions chain-of-custody.",
    },
}


@ChemMCPManager.register_tool
class SamplePreservationGuide(BaseTool):
    """
    样品保存条件推荐：温度、容器、添加剂。
    基于样品类型和目标分析物提供完整的保存方案。
    """
    __version__ = "0.1.0"
    name = "SamplePreservationGuide"
    func_name = "recommend_preservation"
    description = "Recommend optimal sample preservation conditions including temperature, container type, additives, and maximum holding time."
    implementation_description = "Uses a comprehensive database of preservation protocols for different sample types and target analytes based on EPA/ISO guidelines."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Sample Preservation", "Storage", "Quality Assurance", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("sample_type", "str", "N/A", "Type of sample: 'water_inorganic', 'water_organic', 'soil_sediment', 'biological_fluid', 'food_agricultural'."),
        ("target_analytes", "str", "general", "Target analytes or class (e.g., 'metals', 'VOCs', 'pesticides', 'nutrients', 'general')."),
        ("storage_duration_days", "float", "7.0", "Expected storage duration before analysis (days)."),
        ("need_transport", "bool", "False", "Whether samples need to be transported before analysis."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'sample_type target_analytes [storage_duration_days] [need_transport]'"),
    ]

    output_sig = [
        ("container_recommendation", "str", "Recommended container type."),
        ("temperature", "str", "Required storage temperature."),
        ("holding_time_info", "dict", "Maximum holding time for specific analytes."),
        ("preservative_instructions", "list", "Step-by-step preservative addition instructions."),
        ("quality_notes", "list", "QA/QC notes and common pitfalls."),
    ]

    examples = [
        {
            "code_input": {
                "sample_type": "water_inorganic",
                "target_analytes": "metals",
                "storage_duration_days": 30.0,
            },
            "text_input": {
                "input_params": "water_inorganic metals 30.0",
            },
            "output": {
                "container_recommendation": "HDPE bottle (pre-cleaned)",
                "temperature": "4°C (refrigerated)",
                "holding_time_info": {"metals": "6 months with HNO3 pH<2"},
            },
        },
        {
            "code_input": {
                "sample_type": "biological_fluid",
                "target_analytes": "drugs",
                "need_transport": True,
            },
            "text_input": {
                "input_params": "biological_fluid drugs 7.0 True",
            },
            "output": {
                "container_recommendation": "Polypropylene tube",
                "temperature": "-20°C (short-term); -80°C (long-term)",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        sample_type: str,
        target_analytes: str = "general",
        storage_duration_days: float = 7.0,
        need_transport: bool = False,
    ) -> dict:
        """Core logic: recommend preservation conditions."""
        st_key = sample_type.lower().strip()
        analytes_lower = target_analytes.lower().strip()

        # Look up data
        if st_key not in PRESERVATION_DATA:
            matched = None
            for k in PRESERVATION_DATA:
                if st_key in k or k in st_key:
                    matched = k
                    break
            if matched is None:
                raise ChemMCPError(
                    f"Unknown sample type: '{sample_type}'. "
                    f"Available: {', '.join(PRESERVATION_DATA.keys())}"
                )
            st_key = matched

        pdata = PRESERVATION_DATA[st_key]

        # Find relevant preservative info
        relevant_preservatives = []
        if analytes_lower == "general":
            relevant_preservatives = pdata["preservatives"][:3]
        else:
            for p in pdata["preservatives"]:
                if any(a in p["analyte"].lower() for a in analytes_lower.split(",")) or \
                   any(p["analyte"].lower() in a for a in analytes_lower.split(",")):
                    relevant_preservatives.append(p)
            if not relevant_preservatives:
                relevant_preservatives = pdata["preservatives"][:2]

        # Holding time check
        ht_info = self._get_holding_time(pdata["max_holding_time"], analytes_lower, storage_duration_days)
        time_warning = ""
        if ht_info.get("exceeds_max", False):
            time_warning = f"⚠️ WARNING: Requested storage ({storage_duration_days} days) may exceed recommended maximum holding time ({ht_info['max_days']} days)."

        # Transport considerations
        transport_notes = []
        if need_transport:
            transport_notes.append("Use insulated cooler with ice packs (maintain 4°C) or dry ice (for frozen).")
            transport_notes.append("Pack samples upright with absorbent material to contain leaks.")
            transport_notes.append("Include chain-of-custody documentation.")
            if "voc" in analytes_lower or st_key == "water_organic":
                transport_notes.append("Keep VOC samples cool but NOT frozen — avoid pressure buildup in sealed containers.")

        # QA notes
        qa_notes = self._generate_qa_notes(st_key, analytes_lower, storage_duration_days)

        logger.info(f"Preservation guide for {sample_type}, targets={target_analytes}")
        return {
            "sample_type": pdata["name"],
            "container_recommendation": pdata["recommended_container"],
            "temperature": pdata["temperature"],
            "holding_time_info": ht_info,
            "preservative_instructions": relevant_preservatives,
            "transport_requirements": transport_notes if transport_notes else None,
            "time_warning": time_warning,
            "quality_notes": qa_notes,
            "general_notes": pdata["notes"],
        }

    def _get_holding_time(self, ht_dict, analytes, requested_days):
        """Get relevant holding time information."""
        result = {}
        max_allowed = float("inf")

        if analytes == "general":
            result["general"] = ht_dict.get("default", "Check specific analyte class")
            try:
                # Extract numeric value from default string
                import re
                nums = re.findall(r'[\d.]+', str(ht_dict.get("default", "")))
                if nums:
                    max_allowed = min(max_allowed, float(nums[0]) * (365 if "year" in str(ht_dict.get("default", "")).lower() else 1))
                    if "month" in str(ht_dict.get("default", "")).lower():
                        max_allowed = float(nums[0]) * 30
            except (ValueError, IndexError):
                pass
        else:
            for key in ht_dict:
                if key != "default":
                    if any(a in key.lower() for a in analytes.replace(",", " ").split()) or \
                       any(key.lower() in a for a in analytes.replace(",", " ").split()):
                        result[key] = ht_dict[key]
                        try:
                            import re
                            nums = re.findall(r'[\d.]+', str(ht_dict[key]))
                            if nums:
                                val = float(nums[0])
                                if "month" in str(ht_dict[key]).lower():
                                    val *= 30
                                elif "year" in str(ht_dict[key]).lower():
                                    val *= 365
                                max_allowed = min(max_allowed, val)
                        except (ValueError, IndexError):
                            pass

            if not result:
                result["general"] = ht_dict.get("default", "See method-specific guidance")

        result["max_days"] = max_allowed if max_allowed != float("inf") else None
        result["exceeds_max"] = max_allowed < requested_days if max_allowed != float("inf") else False
        return result

    def _generate_qa_notes(self, st_key, analytes, duration):
        notes = [
            "Label all containers with: sample ID, date/time, location, collector name, preservative added.",
            "Use field blanks and trip blanks to monitor contamination during sampling and transport.",
            "Document any deviations from the preservation protocol.",
        ]
        if duration > 30:
            notes.append(f"For long-term storage (>30 days): consider splitting into aliquots to avoid repeated freeze-thaw cycles.")
        if st_key == "water_inorganic":
            notes.append("For metal analysis: use trace-metal grade acids and pre-cleaned containers (10% HNO3 soak).")
        if st_key == "water_organic":
            notes.append("Do NOT use Teflon tape on VOC containers — it can interfere with analysis.")
        if "bio" in st_key:
            notes.append("Process biological samples quickly — enzymatic activity continues even at 4°C.")
        return notes

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            stype = parts[0]
            analytes = parts[1] if len(parts) > 1 else "general"
            days = float(parts[2]) if len(parts) > 2 else 7.0
            transport = parts[3].lower() == "true" if len(parts) > 3 else False
            return self._run_base(stype, analytes, days, transport)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'sample_type analytes [days] [transport_bool]'")
