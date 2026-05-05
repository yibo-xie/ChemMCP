# Log #309: SamplePreservationGuide (样品保存指南)

## Tool Info
- **Tool ID**: 309
- **Class Name**: `SamplePreservationGuide`
- **Module**: `sample_preservation_guide`
- **Version**: 0.1.0

## Core Logic
根据样品类型和目标分析物推荐最佳保存条件（容器、温度、保存时间、添加剂/防腐剂）。内置 EPA/ISO 标准数据库，覆盖水样（无机/有机/VOCs）、生物样品、土壤沉积物等基质的保存规范。

## Input Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| sample_type | str | N/A | 样品类型 (water_inorganic/water_organic/biological_fluid/soil_sediment/tissue/food) |
| target_analytes | str | N/A | 目标分析物 (metals/nutrients/VOCs/SVOCs/drugs/pesticides/cyanide) |
| storage_days | float | 7.0 | 计划保存天数 |
| needs_transport | bool | False | 是否需要运输 |

## Output Fields
| Field | Type | Description |
|-------|------|-------------|
| container_recommendation | str | 推荐容器类型 |
| temperature | str | 保存温度要求 |
| holding_time_info | dict | 最大保存时间（按分析物分类） |
| preservative_instructions | list | 防腐剂添加步骤 |
| transport_requirements | list | 运输要求（如需要） |
| quality_notes | list | QA/QC 注意事项 |

## Example: Water Inorganic / Metals / 7 days
```json
{
  "container_recommendation": "HDPE bottle (pre-cleaned with 10% HNO3)",
  "temperature": "4°C (refrigerated)",
  "holding_time_info": {
    "metals": "6 months (acidified)",
    "Hg": "28 days (acidified, glass container)",
    "Cr(VI)": "24 hours"
  },
  "preservative_instructions": [
    {"step": 1, "additive": "HNO3", "concentration": "pH < 2", "action": "Add concentrated HNO3 to pH < 2"}
  ],
  "quality_notes": ["Use only certified trace-metal grade acid", "Double-cap for transport"]
}
```

## Cherry Studio Config Key: `"--tools", "SamplePreservationGuide"`
