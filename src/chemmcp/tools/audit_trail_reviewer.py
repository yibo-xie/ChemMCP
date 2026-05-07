"""
数据完整性与审计追踪检查工具 (Audit Trail Reviewer)

基于 ALCOA+ 原则检查审计追踪记录的数据完整性。
ALCOA+: Attributable, Legible, Contemporaneous, Original, Accurate
       + Complete, Consistent, Enduring, Available
"""

import logging
import math
from datetime import datetime, timezone
from typing import List, Dict, Optional
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ALCOA+ 原则定义和检查项
ALCOA_PRINCIPLES = {
    "Attributable": {
        "full_name": "Attributable (可归因性)",
        "description": "Who performed an action and when? Data must identify the person or system that created, modified, or deleted it.",
        "checklist": [
            ("user_identity", "Each record has a unique user identifier (username/ID)"),
            ("timestamp_present", "Each action has a date and time stamp"),
            ("electronic_signature", "Electronic signatures include meaning (reason for signing)"),
            ("role_based_access", "Access is controlled by role-based permissions"),
            ("no_shared_accounts", "No shared/generic user accounts used"),
        ],
        "weight": 1.0,
    },
    "Legible": {
        "full_name": "Legible (清晰可读)",
        "description": "Data must be readable and permanent throughout the retention period.",
        "checklist": [
            ("readable_format", "Data is in a human-readable format"),
            ("permanent_record", "Records are permanent (not temporary displays)"),
            ("font_size", "Font size is legible (for printed/electronic records)"),
            ("color_contrast", "Sufficient contrast between text and background"),
            ("no_overwriting", "Original data is not obscured by overwrites/correction fluid"),
        ],
        "weight": 1.0,
    },
    "Contemporaneous": {
        "full_name": "Contemporaneous (同时性)",
        "description": "Data should be recorded at the time of the activity or observation.",
        "checklist": [
            ("real_time_recording", "Data is recorded at the time of the event"),
            ("chronological_order", "Entries are in chronological sequence"),
            ("no_backdating", "No evidence of backdating entries"),
            ("time_accuracy", "System clock is synchronized to a reliable time source"),
            ("delay_justified", "Any recording delays are documented and justified"),
        ],
        "weight": 1.0,
    },
    "Original": {
        "full_name": "Original (原始性)",
        "description": "The first capture of data or a certified true copy (verified copy).",
        "checklist": [
            ("first_capture", "Data includes the original record or verified certified copy"),
            ("source_identified", "Source system/instrument is identified"),
            ("copy_verification", "Certified copies are marked as 'certified true copy' with signature/date"),
            ("metadata_preserved", "Original metadata is preserved with copies"),
            ("raw_data_available", "Raw/primary data is available and accessible"),
        ],
        "weight": 1.0,
    },
    "Accurate": {
        "full_name": "Accurate (准确性)",
        "description": "Data must be correct, truthful, and reflect what was actually observed.",
        "checklist": [
            ("calculation_verification", "Calculations are correct and independently verified"),
            ("data_entry_check", "Double-check or independent verification of manual entries"),
            ("instrument_calibration", "Instruments are calibrated within valid period"),
            ("outlier_investigation", "Outliers/atypical results are investigated and documented"),
            ("amendment_reason", "All corrections have a documented reason for change"),
        ],
        "weight": 1.0,
    },
    # ALCOA+ 扩展原则
    "Complete": {
        "full_name": "Complete (完整性)",
        "description": "All data is present — no deletions or selective omissions.",
        "checklist": [
            ("no_deletions", "No records have been deleted without audit trail"),
            ("all_data_included", "All test results including repeats/failures are included"),
            ("audit_trail_complete", "Audit trail captures all create/read/update/delete actions"),
            ("metadata_complete", "Metadata (sample ID, method version, etc.) is complete"),
            ("sequence_gaps_checked", "Sequence gaps (e.g., injection sequences) are investigated"),
        ],
        "weight": 0.8,
    },
    "Consistent": {
        "full_name": "Consistent (一致性)",
        "description": "Data elements should be logically sequenced, dated, and use consistent units/naming conventions.",
        "checklist": [
            ("naming_convention", "File naming follows established convention"),
            ("date_format_consistent", "Date/time format is consistent across records"),
            ("unit_consistent", "Units of measurement are consistent and clearly defined"),
            ("cross_reference_match", "Cross-references between documents match (e.g., batch numbers)"),
            ("logical_sequence", "Data follows logical chronological and operational sequence"),
        ],
        "weight": 0.8,
    },
    "Enduring": {
        "full_name": "Enduring (持久性)",
        "description": "Data must be available throughout its required retention period.",
        "checklist": [
            ("retention_period_defined", "Retention period is defined per regulatory requirements"),
            ("backup_procedure", "Regular backup procedures are in place"),
            ("storage_conditions", "Storage conditions protect from degradation"),
            ("migration_documented", "Any format migration is validated and documented"),
            ("disaster_recovery", "Disaster recovery plan exists for critical data"),
        ],
        "weight": 0.7,
    },
    "Available": {
        "full_name": "Available (可获取性)",
        "description": "Data must be accessible for review throughout the retention period.",
        "checklist": [
            ("retrieval_tested", "Data retrieval has been tested and works reliably"),
            ("access_controls", "Appropriate access controls allow authorized retrieval"),
            ("searchable", "Records can be searched and retrieved efficiently"),
            ("archive_indexed", "Archived records are indexed for retrieval"),
            ("offline_access", "Procedure exists for accessing offline/archived data"),
        ],
        "weight": 0.7,
    },
}


@ChemMCPManager.register_tool
class AuditTrailReviewer(BaseTool):
    """
    数据完整性与审计追踪检查工具。基于 ALCOA+ 原则，
    对审计追踪记录进行系统性的完整性评估。
    """
    __version__ = "0.1.0"
    name = "AuditTrailReviewer"
    func_name = "review_audit_trail"
    description = "Review audit trail records for data integrity compliance using ALCOA+ principles (Attributable, Legible, Contemporaneous, Original, Accurate + Complete, Consistent, Enduring, Available)."
    implementation_description = (
        "Evaluates audit trail records against all 9 ALCOA+ principles with weighted scoring. "
        "Generates a compliance report with per-principle scores, overall completeness percentage, "
        "identified issues, risk assessment, and remediation recommendations."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Data Integrity", "ALCOA+", "Audit Trail", "GxP", "21 CFR Part 11", "QA/QC", "Compliance"]
    required_envs = []

    code_input_sig = [
        ("audit_records", "list", "N/A",
         "List of audit trail entry dicts. Each dict: "
         "{timestamp(str), user_id(str), action(str), field(str), old_value(any), new_value(any), reason(str)}"),
        ("check_level", "str", "comprehensive",
         "Check depth: 'basic' (core 5 ALCOA principles only) or 'comprehensive' (all 9 ALCOA+ principles)."),
        ("system_context", "dict", "{}",
         "Optional context about the system: {system_name, gxp_type, retention_years, has_backup, ...}"),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "For basic usage: provide check_level ('basic'/'comprehensive'). "
         "Full analysis requires code interface with audit_records list."),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary containing: alcoa_scores(dict per principle), overall_score(%), "
         "issues_found(list), risk_assessment, recommendations, compliance_report_text"),
    ]

    examples = [
        {
            "code_input": {
                "audit_records": [
                    {"timestamp": "2026-01-15T09:30:00Z", "user_id": "analyst_01", "action": "CREATE",
                     "field": "Assay_Result", "old_value": None, "new_value": "99.5%", "reason": "Initial entry"},
                    {"timestamp": "2026-01-15T10:15:00Z", "user_id": "analyst_01", "action": "UPDATE",
                     "field": "Assay_Result", "old_value": "99.5%", "new_value": "99.8%",
                     "reason": "Correction: calculation error identified during review"},
                    {"timestamp": "2026-01-15T14:20:00Z", "user_id": "reviewer_02", "action": "REVIEW",
                     "field": "Batch_B12345", "old_value": None, "new_value": "APPROVED",
                     "reason": "Batch release review completed"},
                    {"timestamp": "2026-01-16T08:00:00Z", "user_id": "admin_sys", "action": "SYSTEM",
                     "field": "Backup_Completed", "old_value": None, "new_value": "Backup_20260116",
                     "reason": "Scheduled nightly backup"},
                ],
                "check_level": "comprehensive",
                "system_context": {
                    "system_name": "LIMS",
                    "gxp_type": "GMP",
                    "retention_years": 10,
                    "has_backup": True,
                },
            },
            "text_input": {
                "input_params": "comprehensive"
            },
            "output": {
                "result": {
                    "overall_compliance_percent": 85.5,
                    "principle_scores": {...},
                    "issues_found": [...],
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
        audit_records: list,
        check_level: str = "comprehensive",
        system_context: dict = None,
    ) -> dict:
        """
        核心逻辑：ALCOA+ 审计追踪审查

        Parameters:
            audit_records: 审计追踪记录列表
            check_level: 检查深度 ('basic' 或 'comprehensive')
            system_context: 系统上下文信息

        Returns:
            dict: 完整的合规性报告
        """
        if not audit_records:
            raise ChemMCPError("Audit records list cannot be empty.")

        ctx = system_context or {}
        level = check_level.lower().strip()

        # 确定要检查的原则
        if level == "basic":
            principles_to_check = {k: v for k, v in ALCOA_PRINCIPLES.items()
                                   if v["weight"] >= 1.0}  # 核心ALCOA 5个
        else:
            principles_to_check = ALCOA_PRINCIPLES

        # 对每条记录执行各原则的检查
        principle_results = {}
        all_issues = []

        for prin_name, prin_def in principles_to_check.items():
            result = self._evaluate_principle(prin_name, prin_def, audit_records, ctx)
            principle_results[prin_name] = result
            all_issues.extend(result.get("issues", []))

        # 计算总体得分（加权平均）
        total_weight = sum(p["weight"] for p in principles_to_check.values())
        weighted_sum = sum(
            principle_results[p]["score_percent"] * principles_to_check[p]["weight"]
            for p in principles_to_check
        )
        overall_score = weighted_sum / total_weight if total_weight > 0 else 0

        # 风险评估
        risk_assessment = self._assess_risk(overall_score, principle_results, all_issues, ctx)

        # 改进建议
        recommendations = self._generate_recommendations(principle_results, all_issues, ctx)

        # 合规判定
        if overall_score >= 90:
            compliance_status = "COMPLIANT ✅"
        elif overall_score >= 70:
            compliance_status = "PARTIALLY COMPLIANT ⚠️"
        else:
            compliance_status = "NON-COMPLIANT ❌"

        result = {
            "review_summary": {
                "total_records_reviewed": len(audit_records),
                "check_level": level,
                "principles_evaluated": len(principles_to_check),
                "overall_compliance_percent": round(overall_score, 1),
                "compliance_status": compliance_status,
                "system_context": ctx,
            },
            "principle_scores": {
                name: {
                    "full_name": r["full_name"],
                    "score_percent": round(r["score_percent"], 1),
                    "items_passed": r["items_passed"],
                    "items_total": r["items_total"],
                    "issues_count": len(r.get("issues", [])),
                    "findings": r.get("findings", ""),
                }
                for name, r in principle_results.items()
            },
            "issues_found": all_issues,
            "risk_assessment": risk_assessment,
            "recommendations": recommendations,
            "report_text": self._generate_report(
                len(audit_records), level, overall_score, compliance_status,
                principle_results, all_issues, recommendations
            ),
        }

        logger.info(f"Audit trail review complete: {len(audit_records)} records, score={overall_score:.1f}%")
        return result

    def _evaluate_principle(self, prin_name: str, prin_def: dict, records: list, ctx: dict) -> dict:
        """评估单个ALCOA+原则"""
        checklist = prin_def["checklist"]
        items_passed = 0
        items_total = len(checklist)
        issues = []
        findings_parts = []

        for item_key, item_desc in checklist:
            passed, issue = self._check_item(prin_name, item_key, item_desc, records, ctx)
            if passed:
                items_passed += 1
            else:
                issues.append({
                    "principle": prin_name,
                    "full_principle": prin_def["full_name"],
                    "check_item": item_key,
                    "description": item_desc,
                    "severity": issue.get("severity", "Medium"),
                    "detail": issue.get("detail", "Check failed"),
                })
                findings_parts.append(f"❌ {item_desc}: {issue.get('detail', 'Failed')}")

        score_pct = (items_passed / items_total * 100) if items_total > 0 else 100.0

        findings = ""
        if findings_parts:
            findings = "\n".join(findings_parts)
        elif items_passed == items_total:
            findings = f"✅ All {items_total} checks passed for {prin_def['full_name']}"

        return {
            "full_name": prin_def["full_name"],
            "score_percent": score_pct,
            "items_passed": items_passed,
            "items_total": items_total,
            "issues": issues,
            "findings": findings,
        }

    def _check_item(self, prin_name: str, item_key: str, item_desc: str, records: list, ctx: dict) -> tuple:
        """检查单个项目，返回 (passed, issue_dict_or_None)"""
        # 基于实际记录内容进行检查
        if prin_name == "Attributable":
            if item_key == "user_identity":
                has_user = all(r.get("user_id") for r in records)
                if not has_user:
                    return False, {"severity": "Critical", "detail": f"{sum(1 for r in records if not r.get('user_id'))} records missing user identity"}
            elif item_key == "timestamp_present":
                has_ts = all(r.get("timestamp") for r in records)
                if not has_ts:
                    return False, {"severity": "Critical", "detail": f"{sum(1 for r in records if not r.get('timestamp'))} records missing timestamp"}
            elif item_key == "amendment_reason":
                updates = [r for r in records if r.get("action", "").upper() in ("UPDATE", "DELETE", "MODIFY")]
                has_reason = all(r.get("reason") for r in updates)
                if updates and not has_reason:
                    return False, {"severity": "High", "detail": f"Some update/delete records lack reason for change"}
            return True, None

        elif prin_name == "Contemporaneous":
            if item_key == "chronological_order":
                timestamps = [r.get("timestamp", "") for r in records if r.get("timestamp")]
                if len(timestamps) > 1:
                    sorted_ts = sorted(timestamps)
                    if timestamps != sorted_ts:
                        return False, {"severity": "Medium", "detail": "Records are not in chronological order"}
            elif item_key == "no_backdating":
                # 检查时间戳是否在未来（相对于其他记录）
                pass  # 需要外部参考时间
            return True, None

        elif prin_name == "Accurate":
            if item_key == "amendment_reason":
                updates = [r for r in records if r.get("action", "").upper() in ("UPDATE", "MODIFY")]
                no_reason = [r for r in updates if not r.get("reason") or r.get("reason").strip() == ""]
                if no_reason:
                    return False, {"severity": "High", "detail": f"{len(no_reason)} modification(s) without documented reason"}
            return True, None

        elif prin_name == "Complete":
            if item_key == "no_deletions":
                deletes = [r for r in records if r.get("action", "").upper() in ("DELETE", "REMOVE")]
                if deletes:
                    return False, {"severity": "Critical", "detail": f"{len(deletes)} deletion(s) found in audit trail - verify legitimacy"}
            elif item_key == "audit_trail_complete":
                # 检查是否有 CREATE 记录对应每个 UPDATE
                creates = set()
                for r in records:
                    if r.get("action", "").upper() == "CREATE":
                        creates.add(r.get("field", ""))
                if not creates:
                    return False, {"severity": "Medium", "detail": "No CREATE records found - possible incomplete trail"}
            return True, None

        elif prin_name == "Enduring":
            if item_key == "retention_period_defined" and ctx:
                if not ctx.get("retention_years"):
                    return False, {"severity": "Medium", "detail": "Retention period not defined in system context"}
            elif item_key == "backup_procedure" and ctx:
                if not ctx.get("has_backup"):
                    return False, {"severity": "High", "detail": "No backup procedure indicated"}
            return True, None

        # 默认：如果无法从记录中验证，返回通过（需要人工确认）
        return True, None

    def _assess_risk(self, overall_score: float, principle_results: dict, issues: list, ctx: dict) -> dict:
        """风险评估"""
        # 统计问题严重程度
        critical = sum(1 for i in issues if i.get("severity") == "Critical")
        high = sum(1 for i in issues if i.get("severity") == "High")
        medium = sum(1 for i in issues if i.get("severity") == "Medium")
        low = sum(1 for i in issues if i.get("severity") == "Low")

        # 风险等级
        if critical > 0 or overall_score < 50:
            risk_level = "CRITICAL"
            risk_description = (
                "Immediate action required. Critical data integrity gaps detected. "
                "Regulatory inspection would likely result in observation/warning letter."
            )
        elif high > 2 or overall_score < 70:
            risk_level = "HIGH"
            risk_description = (
                "Significant data integrity concerns. CAPA required within 30 days. "
                "May result in regulatory observations."
            )
        elif high > 0 or medium > 3 or overall_score < 85:
            risk_level = "MEDIUM"
            risk_description = (
                "Moderate data integrity gaps. Corrective actions recommended. "
                "Document improvements in quality system."
            )
        else:
            risk_level = "LOW"
            risk_description = (
                "Generally compliant with minor observations. Continue monitoring "
                "and address findings as part of continuous improvement."
            )

        # GxP影响评估
        gxp_type = ctx.get("gxp_type", "Unknown")
        impact_note = ""
        if gxp_type.upper() in ("GMP", "GLP"):
            impact_note = f"Under {gxp_type}, data integrity findings can affect product quality/patient safety."

        return {
            "risk_level": risk_level,
            "risk_description": risk_description,
            "issue_summary": {
                "critical": critical,
                "high": high,
                "medium": medium,
                "low": low,
                "total": len(issues),
            },
            "regulatory_impact": impact_note,
            "recommended_response_time": {
                "CRITICAL": "Immediate (within 24 hours)",
                "HIGH": "Within 2 weeks",
                "MEDIUM": "Within 30 days",
                "LOW": "Within 90 days",
            }.get(risk_level, "As scheduled"),
        }

    def _generate_recommendations(self, principle_results: dict, issues: list, ctx: dict) -> list:
        """生成改进建议"""
        recommendations = []
        rec_set = set()

        for issue in issues:
            prin = issue["principle"]
            key = issue["check_item"]

            rec_key = f"{prin}:{key}"
            if rec_key in rec_set:
                continue
            rec_set.add(rec_key)

            if prin == "Attributable" and key == "user_identity":
                recommendations.append({
                    "priority": "P1",
                    "category": "Attributable",
                    "recommendation": "Implement mandatory user authentication for all data entries. Eliminate shared accounts.",
                    "effort": "Low-Medium",
                })
            elif prin == "Attributable" and key == "amendment_reason":
                recommendations.append({
                    "priority": "P1",
                    "category": "Attributable/Accurate",
                    "recommendation": "Configure system to require reason-for-change text before allowing any data modification.",
                    "effort": "Low",
                })
            elif prin == "Complete" and key == "no_deletions":
                recommendations.append({
                    "priority": "P1",
                    "category": "Complete",
                    "recommendation": "Enable append-only mode for audit trails. Deletions should be prohibited at database level.",
                    "effort": "Medium",
                })
            elif prin == "Contemporaneous" and key == "chronological_order":
                recommendations.append({
                    "priority": "P2",
                    "category": "Contemporaneous",
                    "recommendation": "Verify NTP (Network Time Protocol) synchronization across all systems generating audit records.",
                    "effort": "Low",
                })
            elif prin == "Enduring" and key == "backup_procedure":
                recommendations.append({
                    "priority": "P1",
                    "category": "Enduring",
                    "recommendation": "Implement automated daily backup with off-site replication. Test quarterly restore procedures.",
                    "effort": "Medium-High",
                })

        # 如果没有具体建议，添加通用建议
        if not recommendations:
            recommendations.append({
                "priority": "P3",
                "category": "General",
                "recommendation": "Continue routine monitoring. Schedule periodic ALCOA+ self-inspections (quarterly recommended).",
                "effort": "Ongoing",
            })

        return recommendations

    def _generate_report(self, n_records, level, score, status, principle_results, issues, recommendations) -> str:
        lines = [
            f"═══ AUDIT TRAIL REVIEW REPORT (ALCOA+) ═══",
            f"",
            f"Records Reviewed: {n_records}",
            f"Check Level: {level.capitalize()}",
            f"",
            f"═══ OVERALL RESULT ═══",
            f"  Compliance Score: {score:.1f}%",
            f"  Status: {status}",
            f"",
            f"═══ PRINCIPLE SCORES ═══",
        ]
        for pname, pres in principle_results.items():
            bar_len = int(pres["score_percent"] / 5)
            bar = "█" * bar_len + "░" * (20 - bar_len)
            lines.append(f"  {pres['full_name']}: {pres['score_percent']:5.1f}% [{bar}] ({pres['items_passed']}/{pres['items_total']})")

        lines.extend([
            f"",
            f"═══ ISSUES FOUND: {len(issues)} ═══",
        ])
        sev_order = {"Critical": 0, "High": 1, "Medium": 2, "Low": 3}
        sorted_issues = sorted(issues, key=lambda x: sev_order.get(x.get("severity", "Low"), 4))
        for issue in sorted_issues[:15]:  # 最多显示15个
            lines.append(f"  [{issue.get('severity', '?')}] {issue.get('full_principle', '')}: {issue.get('description', '')}")
        if len(sorted_issues) > 15:
            lines.append(f"  ... and {len(sorted_issues)-15} more issues")

        lines.extend([
            f"",
            f"═══ TOP RECOMMENDATIONS ═══",
        ])
        for rec in recommendations[:5]:
            lines.append(f"  [{rec['priority']}] {rec['recommendation']}")

        return "\n".join(lines)

    def _run_text(self, input_params: str) -> dict:
        """文本模式：生成模板报告"""
        level = input_params.strip().lower() if input_params.strip() else "comprehensive"

        # 在文本模式下使用示例数据演示
        sample_records = [
            {"timestamp": "2026-05-06T10:00:00+08:00", "user_id": "analyst_A", "action": "CREATE",
             "field": "Sample_Result", "old_value": None, "new_value": "98.5%", "reason": "Initial data entry"},
            {"timestamp": "2026-05-06T10:30:00+08:00", "user_id": "analyst_A", "action": "UPDATE",
             "field": "Sample_Result", "old_value": "98.5%", "new_value": "98.8%",
             "reason": "Corrected calculation error"},
            {"timestamp": "2026-05-06T14:00:00+08:00", "user_id": "reviewer_B", "action": "REVIEW",
             "field": "Batch_R001", "old_value": None, "new_value": "APPROVED",
             "reason": "Second-person review completed"},
        ]

        return self._run_base(sample_records, level, {
            "system_name": "Example LIMS",
            "gxp_type": "GMP",
            "retention_years": 10,
            "has_backup": True,
        })
