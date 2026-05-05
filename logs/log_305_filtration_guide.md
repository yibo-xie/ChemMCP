# Log #305: FiltrationGuide (滤膜选择指南)

## Tool Info
- **Tool ID**: 305
- **Class Name**: `FiltrationGuide`
- **Module**: `filtration_guide`
- **Version**: 0.1.0

## Core Logic
滤膜选择指导，根据应用场景、孔径、溶剂类型和目标分析物推荐最佳滤膜材质。内置 12 种常见滤膜数据库（CA/Nylon/PES/PTFE/PVDF/RC/MCE/PP/GF/PES-sterile/PTFE-hydrophilic/Nylon-66），包含 pH 范围、耐温性、结合特性等参数，通过评分算法排序推荐。

## Input Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| application | str | N/A | 应用场景 (HPLC_prep/sterile_filtration/particle_analysis/gas_parging/cell_culture) |
| pore_size_um | float | 0.45 | 目标孔径 (μm) |
| solvent_type | str | "aqueous" | 溶剂类型 (aqueous/organic/strong_acid/strong_base/buffer) |
| analyte_type | str | "general" | 分析物类型 (general/protein/peptide/cells) |
| minimize_binding | bool | False | 是否优先低结合 |
| max_temperature_c | float | 25.0 | 最高操作温度 (°C) |

## Output Fields
| Field | Type | Description |
|-------|------|-------------|
| primary_recommendation | dict | 首选滤膜详情（材质、孔径、pH范围、兼容性） |
| alternatives | list | 备选滤膜列表 |
| selection_rationale | str | 选择理由 |
| usage_tips | list | 使用注意事项 |

## Example: HPLC Sample Prep / 0.45μm / Organic → Nylon
```json
{
  "primary_recommendation": {
    "material": "Nylon (Polyamide, NY)",
    "available_pore_sizes_um": [0.1, 0.2, 0.22, 0.45, 0.65, 0.8, 1.0, 2.0, 3.0],
    "recommended_pore_um": 0.45,
    "ph_range": "3-12",
    "max_temperature_c": 180,
    "binding_characteristics": "Moderate protein binding"
  },
  "alternatives": ["PES", "PTFE", "PVDF"]
}
```

## Cherry Studio Config Key: `"--tools", "FiltrationGuide"`
