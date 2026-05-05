# Log #308: MatrixMatchingAdvisor (基质匹配建议)

## Tool Info
- **Tool ID**: 308
- **Class Name**: `MatrixMatchingAdvisor`
- **Module**: `matrix_matching_advisor`
- **Version**: 0.1.0

## Core Logic
根据样品基质和检测技术评估基质效应严重程度，推荐最佳基质匹配策略。内置 8 种常见基质的基质效应数据库（food/water/blood_plasma/urine/soil/tissue/environmental/pharmaceutical），涵盖离子抑制/增强、内标选择、校准策略等。

## Input Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| sample_matrix | str | N/A | 样品基质 (food/water_environmental/blood_plasma/urine/soil/tissue) |
| detection_technique | str | N/A | 检测技术 (LC-MS/GC-MS/HPLC-UV/ICP-MS/ICP-OES) |
| quantification_level | str | "trace_level" | 定量水平 (trace_level/major_component/screening) |
| target_analytes | str | "general" | 目标分析物类型 |

## Output Fields
| Field | Type | Description |
|-------|------|-------------|
| matrix_effect_assessment | dict | 基质效应评估（严重程度、主要来源、影响机制） |
| recommended_strategy | dict | 推荐匹配策略（方法名、原理、操作步骤） |
| internal_standard_advice | dict | 内标选择建议 |
| calibration_strategy | dict | 校准曲线策略 |
| workflow | list | 完整工作流程步骤 |

## Example: Food / LC-MS/MS / Trace Level
```json
{
  "matrix_effect_assessment": {
    "severity": "Very High",
    "primary_source": "Co-eluting matrix components (lipids, pigments, sugars)",
    "mechanism": "Ion suppression in ESI source"
  },
  "recommended_strategy": {
    "method_name": "SPE after extraction",
    "principle": "Remove co-extractives before analysis",
    "steps": ["Extract with appropriate solvent", "SPE cleanup", "Reconstitute in mobile phase"]
  },
  "internal_standard_advice": {
    "recommendation": "SIL-IS essential for accurate quantification",
    "type": "Stable isotope-labeled internal standard"
  },
  "workflow": [
    {"step": 1, "name": "Sample extraction", "action": "..."},
    ...
  ]
}
```

## Cherry Studio Config Key: `"--tools", "MatrixMatchingAdvisor"`
