# Log #304: DigestionProtocolSelector (消解方法选择)

## Tool Info
- **Tool ID**: 304
- **Class Name**: `DigestionProtocolSelector`
- **Module**: `digestion_protocol_selector`
- **Version**: 0.1.0

## Core Logic
根据样品基质和目标元素推荐最佳消解方法（湿法消解、微波消解、干灰化法）。内置元素挥发性数据库，对含 Hg、As 等易挥发元素的样品自动给出警告和特殊处理建议。

## Input Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| sample_type | str | N/A | 样品类型 (food/soil/plant/water/biological/sediment) |
| target_elements | str | N/A | 目标元素，逗号分隔 (如 Pb,Cd,As,Hg) |
| equipment_available | str | "full" | 设备条件 (full/basic/minimal) |
| priority | str | "recovery" | 优先级 (recovery/speed/simplicity) |

## Output Fields
| Field | Type | Description |
|-------|------|-------------|
| recommended_method | dict | 推荐方法详情（名称、中文名、原理） |
| protocol_steps | list | 分步操作流程 |
| warnings | list | 安全警告 |
| alternatives | list | 备选方法 |
| rationale | str | 选择理由 |

## Example: Food / Pb,Cd,As,Hg / Full Equipment → Microwave Digestion
```json
{
  "recommended_method": {
    "name": "Microwave-Assisted Digestion",
    "name_cn": "微波消解",
    "principle": "Closed-vessel high-pressure acid digestion"
  },
  "warnings": ["WARNING: Hg and As are volatile - use reflux condenser or closed vessel"],
  "alternatives": ["Wet Acid Digestion (湿法消解)", "Dry Ashing (干灰化法)"]
}
```

## Cherry Studio Config Key: `"--tools", "DigestionProtocolSelector"`
