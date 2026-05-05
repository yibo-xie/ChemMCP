# Log #310: HomogenizationProtocol (均质化方案)

## Tool Info
- **Tool ID**: 310
- **Class Name**: `HomogenizationProtocol`
- **Module**: `homogenization_protocol`
- **Version**: 0.1.0

## Core Logic
固体样品均质化方案设计，根据样品类型、状态和分析物稳定性推荐最佳均质方法。内置 8 种均质技术数据库（研钵/球磨机/高速匀浆器/珠磨/组织研磨仪/冷冻研磨/超声/旋风磨），包含操作参数、适用场景和关键控制点。

## Input Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| sample_type | str | N/A | 样品类型 (animal_tissue/plant_tissue/soil/food/microbial_cells/pharmaceutical/sediment/environmental) |
| sample_state | str | "fresh" | 样品状态 (fresh/frozen/dried) |
| analyte_stability | str | "stable" | 分析物稳定性 (stable/thermolabile/light_sensitive/oxidation_sensitive) |
| sample_amount_g | float | 10.0 | 样品量 (g) |

## Output Fields
| Field | Type | Description |
|-------|------|-------------|
| recommended_method | dict | 推荐方法详情（名称、设备、原理） |
| protocol_steps | list | 分步操作流程 |
| rationale | str | 选择理由 |
| alternatives | list | 备选方案 |
| critical_points | list | 关键控制点 |

## Example: Animal Tissue / Frozen / Thermolabile → Cryogenic Grinding
```json
{
  "recommended_method": {
    "method": "Mortar & Pestle (Cryogenic Grinding)",
    "equipment": "Liquid N2-cooled mortar and pestle or cryogenic mill",
    "principle": "Sample embrittlement at liquid N2 temperature, mechanical grinding"
  },
  "protocol_steps": [
    {"step": 1, "name": "Pre-cooling", "action": "Cool mortar and pestle with liquid N2"},
    {"step": 2, "name": "Sample addition", "action": "Add frozen tissue (~10g) to cooled mortar"},
    ...
  ],
  "critical_points": ["Keep sample frozen during grinding", "Avoid thawing → analyte degradation", ...]
}
```

## Cherry Studio Config Key: `"--tools", "HomogenizationProtocol"`
