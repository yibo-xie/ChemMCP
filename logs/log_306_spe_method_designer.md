# Log #306: SPEMethodDesigner (固相萃取方法设计)

## Tool Info
- **Tool ID**: 306
- **Class Name**: `SPEMethodDesigner`
- **Module**: `spe_method_designer`
- **Version**: 0.1.0

## Core Logic
固相萃取方法设计，根据分析物性质、样品基质和检测方法自动选择最佳 SPE 填料并生成完整 6 步操作方案（活化→平衡→上样→淋洗→干燥→洗脱）。内置 10 种填料数据库（C18/C8/HLB/WAX/WCX/SAX/SCX/Florisil/PS-DVB/GCB），包含推荐溶剂和体积。

## Input Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| analytes | str | N/A | 目标分析物（逗号分隔） |
| analyte_properties | str | "nonpolar" | 分析物性质 (nonpolar/acidic/basic/polar/mixed) |
| sample_matrix | str | "water" | 样品基质 |
| sample_volume_ml | float | 100.0 | 样品体积 (mL) |
| detection_method | str | "LC-MS" | 检测方法 |
| ph_adjustment | float | None | pH 调节值 |

## Output Fields
| Field | Type | Description |
|-------|------|-------------|
| recommended_sorbent | str | 推荐填料名称及描述 |
| method_protocol | dict | 完整操作方案（cartridge, steps, post_elution） |
| ph_recommendation | dict | pH 调节建议 |
| alternatives | list | 备选填料列表 |
| tips | list | 关键注意事项 |

## Example: Pesticides/Drugs / Water / 500mL / LC-MS → HLB
```json
{
  "recommended_sorbent": "Hydrophilic-Lipophilic Balance (N-vinylpyrrolidone-divinylbenzene copolymer)",
  "method_protocol": {
    "cartridge_recommendation": {"size": "6 mL (200 mg)", "reason": "..."},
    "steps": [
      {"step": 1, "name": "Conditioning (活化)", "action": "Pass 6 mL of Methanol...", "solvent": "Methanol", "volume_ml": 6.0},
      {"step": 2, "name": "Equilibration (平衡)", ...},
      {"step": 3, "name": "Sample Loading (上样)", ...},
      {"step": 4, "name": "Washing (淋洗)", ...},
      {"step": 5, "name": "Drying (干燥)", ...},
      {"step": 6, "name": "Elution (洗脱)", ...}
    ]
  },
  "tips": [
    "Never let sorbent dry completely between conditioning and loading",
    "Flow rates: loading ≤5 mL/min, elution ≤1 mL/min"
  ]
}
```

## Cherry Studio Config Key: `"--tools", "SPEMethodDesigner"`
