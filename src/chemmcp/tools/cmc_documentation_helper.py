"""
CMC 文档编写辅助工具 (CMC Documentation Helper)

提供 Chemistry, Manufacturing, and Controls (CMC) 文档编写指导和模板。
涵盖 ICH M4Q (CTD Module 3 Quality) 各章节。
"""

import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# CTD Module 3 (Quality) - CMC 章节结构（基于 ICH M4Q）
CTD_MODULE_3_STRUCTURE = {
    "S.1_General_Information": {
        "section_code": "S.1",
        "title": "General Information",
        "subsections": {
            "S.1.1": {"name": "Nomenclature", "content_guidance": [
                "Chemical name(s): systematic (IUPAC), common/trade names",
                "CAS Registry Number",
                "INN (International Nonproprietary Name) if applicable",
                "Molecular formula and weight",
                "Structural representation (2D chemical structure)",
            ]},
            "S.1.2": {"name": "Structure", "content_guidance": [
                "Structural formula including relative stereochemistry",
                "Molecular formula: CxHyOzNw...",
                "Molecular weight (average and monoisotopic if relevant)",
            ]},
            "S.1.3": {"name": "General Properties", "content_guidance": [
                "Physical form/appearance (solid, liquid, crystalline, amorphous)",
                "Physicochemical properties: solubility profile (pH-dependent), pKa, partition coefficient (log P)",
                "Ionization constant(s)",
                "Solid-state properties: polymorphism, hygroscopicity",
            ]},
        },
    },
    "S.2_Manufacture": {
        "section_code": "S.2",
        "title": "Manufacture",
        "subsections": {
            "S.2.1": {"name": "Manufacturer(s)", "content_guidance": [
                "Name, address, and responsibility of each manufacturer/site",
                "Contract manufacturers with their specific operations",
                "Flow diagram showing manufacturing sites and supply chain",
            ]},
            "S.2.2": {"name": "Description of Manufacturing Process", "content_guidance": [
                "Process flow diagram (block flow or detailed process flow)",
                "Step-by-step description of the synthetic route / production process",
                "In-process controls at each critical step",
                "Equipment list for key unit operations",
                "Yield expectations at each step",
            ]},
            "S.2.3": {"name": "Control of Materials", "content_guidance": [
                "List of all raw materials, solvents, reagents, catalysts",
                "Specifications and testing for starting materials (SMs)",
                "Justification for SM specifications (especially for critical SMs)",
                "Control of materials with potential carry-over concerns",
            ]},
            "S.2.4": {"name": "Controls of Critical Steps and Intermediates", "content_guidance": [
                "Identification of critical process parameters (CPPs)",
                "In-process controls (IPCs) with acceptance criteria",
                "Intermediate specifications and hold times",
                "Process validation approach summary",
            ]},
            "S.2.5": {"name": "Process Validation", "content_guidance": [
                "Process development summary (lab → pilot → commercial scale)",
                "Validation protocol and report references",
                "CPP ranges established from development studies",
                "Continued Process Verification (CPV) plan",
            ]},
            "S.2.6": {"name": "Manufacturing Process Development", "content_guidance": [
                "Evolution of the manufacturing process from development to commercial",
                "Key experiments that defined process parameters",
                "Risk assessment for process parameters (QbD approach if used)",
                "Scale-up considerations and justification",
            ]},
        },
    },
    "S.3_Characterization": {
        "section_code": "S.3",
        "title": "Characterization",
        "subsections": {
            "S.3.1": {"name": "Elucidation of Structure and Other Characteristics", "content_guidance": [
                "Structural elucidation data: NMR (¹H, ¹³C, 2D), MS, IR, UV-Vis",
                "Elemental analysis results",
                "Stereochemistry determination (chiral HPLC, X-ray, optical rotation)",
                "Solid-state characterization: XRPD, DSC, TGA, SEM (if applicable)",
                "Discussion of polymorphic forms identified",
            ]},
            "S.3.2": {"name": "Impurities", "content_guidance": [
                "Organic impurities: process-related (starting materials, by-products, intermediates, degradants)",
                "Organic impurities: identification and qualification thresholds per ICH Q3A/Q3B",
                "Inorganic impurities: catalysts, heavy metals (ICH Q3D)",
                "Residual solvents: classification and limits per ICH Q3C",
                "Impurity fate mapping (where do they go in the process?)",
                "Justification for impurity specifications (toxicological data if above identification threshold)",
            ]},
        },
    },
    "S.4_Control_of_Drug_Substance": {
        "section_code": "S.4",
        "title": "Control of Drug Substance",
        "subsections": {
            "S.4.1": {"name": "Specification", "content_guidance": [
                "Complete specification table with tests, methods, and acceptance criteria",
                "Justification for each specification parameter",
                "Reference to pharmacopeial monograph if applicable (USP/EP/JP)",
            ]},
            "S.4.2": {"name": "Analytical Procedures", "content_guidance": [
                "Validated analytical method descriptions per ICH Q2(R1)",
                "Method validation summaries (or full reports in appendices)",
                "System suitability requirements",
                "Reference standard characterization (primary vs working standard)",
            ]},
            "S.4.3": {"name": "Validation of Analytical Procedures", "content_guidance": [
                "Summary of validation data for each compendial/non-compendial method",
                "Cross-reference to full validation reports",
                "Method transfer documentation (if applicable)",
            ]},
            "S.4.4": {"name": "Batch Analysis", "content_guidance": [
                "Data from at least 3 pilot/commercial scale batches",
                "Comparison against proposed specification",
                "Discussion of any out-of-specification trends",
            ]},
            "S.4.5": {"name": "Justification of Specification", "content_guidance": [
                "Scientific rationale for each test and acceptance criterion",
                "Batch history data supporting criteria",
                "Regulatory precedent and pharmacopeial alignment",
                "Safety-based justification for impurity limits",
            ]},
        },
    },
    "S.5_Reference_Standards_or_Materials": {
        "section_code": "S.5",
        "title": "Reference Standards or Materials",
        "content_guidance": [
            "Primary reference standard: source, characterization, certificate of analysis",
            "Working reference standard: preparation, assignment of potency, requalification schedule",
            "Handling and storage conditions",
            "Stability of reference standards",
        ],
    },
    "S.6_Container_Closure_System": {
        "section_code": "S.6",
        "title": "Container Closure System",
        "content_guidance": [
            "Description of container/closure system (composition, dimensions)",
            "Qualification data (extractables/leachables if applicable)",
            "Protection provided (light, moisture, oxygen barrier)",
            "Compatibility with drug substance/product",
            "Regulatory status (DMF filed? Pharmacopeial chapter?)",
        ],
    },
    "S.7_Stability": {
        "section_code": "S.7",
        "title": "Stability",
        "subsections": {
            "S.7.1": {"name": "Stability Summary and Conclusions", "content_guidance": [
                "Summary of stability data from long-term, accelerated, and stress conditions",
                "Proposed shelf life and storage conditions with statistical basis",
                "Post-approval stability commitment (commitment protocol)",
            ]},
            "S.7.2": {"name": "Post-Approval Stability Protocol", "content_guidance": [
                "Stability study design for annual batches (per ICH Q1A(R2))",
                "Stability-indicating methods used",
                "Data evaluation approach (trend analysis, statistical modeling)",
            ]},
            "S.7.3": {"name": "Stability Data", "content_guidance": [
                "Tabulated stability data (appearance, assay, degradants, dissolution, etc.)",
                "Graphical presentation of trend data",
                "Statistical analysis (regression, shelf life estimation)",
            ]},
        },
    },
}

# 制剂 (Drug Product) 章节
P2_FORMULATION = {
    "P.2.1": {"name": "Description of Composition", "components": [
        ("Active Pharmaceutical Ingredient (API)", "Name, strength/amount per unit"),
        ("Excipients", "List with function and amount per unit"),
        ("Qualitative/Quantitative composition table", "Full formulation breakdown"),
        ("Overage justification", "If applicable, scientific rationale for excess API"),
    ]},
    "P.2.2": {"name": "Pharmaceutical Development", "components": [
        ("Component selection", "Rationale for excipient choice (functionality, compatibility)"),
        ("Formulation optimization", "DOE studies, critical formulation attributes"),
        ("Drug-excipient compatibility", "Binary/ternary compatibility study results"),
        ("Bioavailability/performance justification", "For modified release products"),
        ("Manufacturability", "Process-formulation interplay discussion"),
    ]},
    "P.2.3": {"name": "Manufacture", "components": [
        ("Unit operation description", "Blending, granulation, compression, coating, filling, etc."),
        ("Equipment", "Type and capacity of key equipment"),
        ("In-process controls", "Critical quality attributes monitored during manufacture"),
    ]},
}


@ChemMCPManager.register_tool
class CMCDocumentationHelper(BaseTool):
    """
    CMC 文档编写辅助工具。提供 ICH M4Q (CTD Module 3 Quality) 格式的
    CMC 文档模板和编写指导，涵盖原料药(DS)和制剂(DP)各章节。
    """
    __version__ = "0.1.0"
    name = "CMCDocumentationHelper"
    func_name = "assist_cmc_documentation"
    description = "Provide structured templates and content guidance for CMC (Chemistry, Manufacturing, Controls) documentation per ICH M4Q (CTD Module 3)."
    implementation_description = (
        "Contains complete CTD Module 3 (Quality) section outlines for Drug Substance (S.1-S.7) "
        "and Drug Product (P.1-P.5) with content guidance, regulatory references, "
        "and data requirements per ICH M4Q(R2), Q6A/B, Q8-Q12 guidelines."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["CMC", "Regulatory Documentation", "CTD", "ICH M4Q", "Pharmaceutical Development", "QA/QC"]
    required_envs = []

    code_input_sig = [
        ("document_section", "str", "all",
         "Specific section code (e.g., 'S.1', 'S.2', 'S.3', 'S.4', 'P.2', etc.) or 'all' for full outline."),
        ("product_type", "str", "DS",
         "'DS' for Drug Substance sections (S.1-S.7) or 'DP' for Drug Product sections (P.1-P.5)."),
        ("drug_name", "str", "",
         "Name of the drug substance or product (used to customize template headers)."),
        ("include_regulatory_refs", "bool", "True",
         "Whether to include ICH guideline references for each section."),
        ("detail_level", "str", "standard",
         "Detail level: 'outline' (section titles only), 'standard' (with subsections), "
         "'detailed' (with full content guidance items)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "String format: '[document_section] [product_type] [drug_name] [detail_level]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary containing: document_structure(dict with sections/subsections), "
         "content_guidance_per_section, regulatory_references, writing_checklist, word_count_estimate"),
    ]

    examples = [
        {
            "code_input": {
                "document_section": "S.4",
                "product_type": "DS",
                "drug_name": "ExampleAPI",
                "include_regulatory_refs": True,
                "detail_level": "detailed",
            },
            "text_input": {
                "input_params": "S.4 DS ExampleAPI detailed"
            },
            "output": {
                "result": {
                    "section_code": "S.4",
                    "section_title": "Control of Drug Substance",
                    "template": {...},
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
        document_section: str = "all",
        product_type: str = "DS",
        drug_name: str = "",
        include_regulatory_refs: bool = True,
        detail_level: str = "standard",
    ) -> dict:
        """
        核心逻辑：生成CMC文档模板

        Parameters:
            document_section: 章节代码 ('all' 或具体章节如 'S.1', 'S.4', 'P.2')
            product_type: 产品类型 ('DS' 或 'DP')
            drug_name: 药物名称
            include_regulatory_refs: 是否包含法规引用
            detail_level: 详细程度

        Returns:
            dict: CMC文档模板和指导
        """
        pt = product_type.upper().strip()
        if pt not in ("DS", "DP"):
            raise ChemMCPError(f"product_type must be 'DS' or 'DP', got '{product_type}'.")

        # 法规引用映射
        reg_refs = self._get_regulatory_references(pt)

        # 根据章节选择内容
        section = document_section.strip()

        if pt == "DS":
            if section.upper() == "ALL":
                selected_sections = CTD_MODULE_3_STRUCTURE
            elif section in CTD_MODULE_3_STRUCTURE:
                selected_sections = {section: CTD_MODULE_3_STRUCTURE[section]}
            else:
                # 模糊匹配
                matched = None
                for key in CTD_MODULE_3_STRUCTURE:
                    if section.replace(".", "_") in key or key.startswith(section):
                        matched = key
                        break
                if matched:
                    selected_sections = {matched: CTD_MODULE_3_STRUCTURE[matched]}
                else:
                    available = ", ".join(CTD_MODULE_3_STRUCTURE.keys())
                    raise ChemMCPError(f"Unknown section '{section}'. Available DS sections: {available}")
        else:  # DP
            selected_sections = self._get_dp_sections(section)

        # 构建输出文档结构
        doc_structure = {}
        total_items = 0
        writing_checklist = []

        for sec_key, sec_data in selected_sections.items():
            sec_entry = self._process_section(
                sec_key, sec_data, drug_name, detail_level,
                include_regulatory_refs, reg_refs, writing_checklist
            )
            doc_structure[sec_key] = sec_entry
            total_items += sec_entry.get("item_count", 0)

        result = {
            "document_info": {
                "type": f"{'Drug Substance' if pt == 'DS' else 'Drug Product'} CMC Documentation",
                "framework": "ICH M4Q(R2) - Common Technical Document Module 3 (Quality)",
                "section_requested": section,
                "sections_included": list(selected_sections.keys()),
                "drug_name": drug_name or "[DRUG NAME]",
                "detail_level": detail_level,
            },
            "document_structure": doc_structure,
            "regulatory_references": reg_refs if include_regulatory_refs else {},
            "writing_checklist": writing_checklist,
            "total_content_items": total_items,
            "report_text": self._generate_report(pt, section, drug_name, total_items, doc_structure),
        }

        logger.info(f"Generated CMC documentation template: {pt}, section={section}, {total_items} items")
        return result

    def _get_regulatory_references(self, pt: str) -> dict:
        """获取相关法规引用"""
        refs = {
            "Core Guidelines": {
                "ICH M4Q(R2)": "Common Technical Document — Module 3: Quality",
                "ICH Q6A": "Specifications: Test Procedures and Acceptance Criteria for New Drug Substances and Products: Chemical Substances",
                "ICH Q6B": "Specifications: Test Procedures and Acceptance Criteria for Biotechnological/Biological Products",
            },
            "Quality Guidelines": {
                "ICH Q2(R1)": "Validation of Analytical Procedures",
                "ICH Q3A(R2)": "Impurities in New Drug Substances",
                "ICH Q3B(R2)": "Impurities in New Drug Products",
                "ICH Q3C(R2)": "Guideline on Residual Solvents",
                "ICH Q3D": "Guideline for Elemental Impurities",
                "ICH Q6A/Q6B": "Specifications",
                "ICH Q7": "GMP for Active Pharmaceutical Ingredients",
                "ICH Q8(R2)": "Pharmaceutical Development",
                "ICH Q9": "Quality Risk Management",
                "ICH Q10": "Pharmaceutical Quality System",
                "ICH Q11": "Development and Manufacture of Drug Substances",
                "ICH Q12": "Lifecycle Management",
            },
            "Stability Guidelines": {
                "ICH Q1A(R2)": "Stability Testing of New Drug Substances and Products",
                "ICH Q1B": "Photostability Testing",
                "ICH Q1C": "Stability Testing for New Dosage Forms",
                "ICH Q1E": "Stability Data Evaluation",
            },
        }
        return refs

    def _get_dp_sections(self, section: str) -> dict:
        """获取制剂章节"""
        dp_structure = {
            "P.1_Description": {
                "section_code": "P.1",
                "title": "Description and Composition of Drug Product",
                "content_guidance": [
                    "Product description (dosage form, strength, route of administration)",
                    "List of components (quantitative and qualitative composition)",
                    "Pharmacopeial compliance statement",
                ],
            },
            "P.2_Pharmaceutical_Development": {
                "section_code": "P.2",
                "title": "Pharmaceutical Development",
                "content_guidance": [
                    "Drug product formulation composition",
                    "Manufacturing process development",
                    "Container closure system",
                    "Compatibility studies",
                    "Critical quality attributes (CQAs)",
                ],
            },
            "P.3_Manufacture": {
                "section_code": "P.3",
                "title": "Manufacture",
                "content_guidance": [
                    "Manufacturers",
                    "Batch formula",
                    "Description of manufacturing process",
                    "Controls of in-process materials",
                    "Process validation",
                ],
            },
            "P.4_Control_of_Excipients": {
                "section_code": "P.4",
                "title": "Control of Excipients",
                "content_guidance": [
                    "Excipient specifications",
                    "Analytical procedures",
                    "Justification",
                ],
            },
            "P.5_Control_of_Drug_Product": {
                "section_code": "P.5",
                "title": "Control of Drug Product",
                "content_guidance": [
                    "Specification (tests, methods, acceptance criteria)",
                    "Analytical procedure validation",
                    "Batch analysis data",
                    "Justification of specification",
                ],
            },
        }

        if section.upper() == "ALL":
            return dp_structure
        elif section in dp_structure:
            return {section: dp_structure[section]}
        else:
            available = ", ".join(dp_structure.keys())
            raise ChemMCPError(f"Unknown DP section '{section}'. Available: {available}")

    def _process_section(self, sec_key, sec_data, drug_name, detail_level,
                         include_refs, reg_refs, checklist) -> dict:
        """处理单个章节"""
        item_count = 0
        processed_subsections = {}

        if "subsections" in sec_data:
            for sub_key, sub_data in sec_data["subsections"].items():
                if detail_level == "outline":
                    sub_output = {"name": sub_data["name"]}
                elif detail_level == "detailed":
                    sub_output = {
                        "name": sub_data["name"],
                        "content_guidance": sub_data.get("content_guidance", []),
                        "regulatory_reference": self._find_relevant_ref(sub_key, reg_refs) if include_refs else None,
                    }
                    item_count += len(sub_data.get("content_guidance", []))
                else:  # standard
                    sub_output = {
                        "name": sub_data["name"],
                        "num_content_items": len(sub_data.get("content_guidance", [])),
                    }
                processed_subsections[sub_key] = sub_output
                checklist.append(f"☐ Write {sec_data.get('section_code', sec_key)}.{sub_key}: {sub_data['name']}")

        elif "content_guidance" in sec_data:
            if detail_level != "outline":
                item_count += len(sec_data["content_guidance"])

        return {
            "section_code": sec_data.get("section_code", sec_key),
            "title": sec_data.get("title", sec_key),
            "subsections": processed_subsections if processed_subsections else None,
            "content_guidance": sec_data.get("content_guidance") if detail_level == "detailed" else None,
            "item_count": item_count,
        }

    def _find_relevant_ref(self, sub_key: str, refs: dict) -> str:
        """查找相关的法规引用"""
        ref_mapping = {
            "S.4.1": "ICH Q6A/Q6B",
            "S.4.2": "ICH Q2(R1)",
            "S.4.3": "ICH Q2(R1)",
            "S.3.2": "ICH Q3A/Q3B/Q3C/Q3D",
            "S.7.1": "ICH Q1A(R2)/Q1E",
            "S.7.2": "ICH Q1A(R2)",
            "S.7.3": "ICH Q1A(R2)/Q1E",
            "S.2.2": "ICH Q11/Q8",
            "S.2.5": "ICH Q7",
        }
        return ref_mapping.get(sub_key, "See ICH M4Q(R2)")

    def _generate_report(self, pt, section, drug_name, total_items, doc_structure) -> str:
        lines = [
            f"═══ CMC DOCUMENTATION TEMPLATE ═══",
            f"",
            f"Document Type: {'Drug Substance (DS)' if pt == 'DS' else 'Drug Product (DP)'}",
            f"Framework: ICH M4Q(R2) - CTD Module 3 (Quality)",
            f"Drug Name: {drug_name or '[DRAG NAME]'}",
            f"Section: {section}",
            f"Total Content Items: {total_items}",
            f"",
            f"─── Sections Included ───",
        ]
        for key, val in doc_structure.items():
            title = val.get("title", key)
            subs = val.get("subsections", {})
            if subs:
                lines.append(f"\n  {key}: {title}")
                for sk, sv in subs.items():
                    sname = sv.get("name", "") if isinstance(sv, dict) else str(sv)
                    lines.append(f"    └─ {sk}: {sname}")
            else:
                lines.append(f"  {key}: {title}")

        lines.append(f"\nUse writing_checklist in output for step-by-step tracking.")
        return "\n".join(lines)

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            section = parts[0].upper() if parts else "ALL"
            pt = parts[1].upper() if len(parts) > 1 else "DS"
            name = parts[2] if len(parts) > 2 else ""
            level = parts[3] if len(parts) > 3 else "standard"

            return self._run_base(section, pt, name, True, level)
        except (ValueError, IndexError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
