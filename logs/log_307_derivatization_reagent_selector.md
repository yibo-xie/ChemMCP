# Log #307: DerivatizationReagentSelector (衍生化试剂选择)

## Tool Info
- **Tool ID**: 307
- **Class Name**: `DerivatizationReagentSelector`
- **Module**: `derivatization_reagent_selector`
- **Version**: 0.1.0

## Core Logic
根据官能团类型和检测方法选择最佳衍生化试剂。内置 12 种常用衍生化试剂数据库（BSTFA/MSTFA/TMS/PFBBr/DNPH/FMOC-Cl/dansyl-Cl/BF3-MeOH/TMO/AA/PAA/AC），包含反应条件（温度、时间、溶剂）、优缺点和适用检测方法。

## Input Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| functional_group | str | N/A | 目标官能团 (-OH, -NH2, -COOH, -SH, C=O, -COOH) |
| detection_method | str | N/A | 检测方法 (GC-MS, GC-FID, HPLC-UV, HPLC-FLD, LC-MS) |
| constraint | str | "general" | 约束条件 (thermal_stability/safety/speed/sensitivity) |

## Output Fields
| Field | Type | Description |
|-------|------|-------------|
| recommended_reagent | dict | 推荐试剂详情（名称、反应类型、温度、时间、溶剂、机理、优缺点） |
| protocol | list | 分步操作流程 |
| alternatives | list | 备选试剂列表 |

## Example: -COOH / GC-FID / Thermal Stability → BSTFA+TMCS or TMO
```json
{
  "recommended_reagent": {
    "reagent": "BSTFA + TMCS (1%)",
    "full_name": "N,O-Bis(trimethylsilyl)trifluoroacetamide + Trimethylchlorosilane",
    "reaction_type": "Silylation (TMS ether/ester)",
    "reaction_temp_c": "60-80",
    "reaction_time_min": "15-30",
    "solvent": "Pyridine, acetonitrile, DMF",
    "pros": ["Single-step", "Mild conditions", "By-products volatile"],
    "cons": ["Moisture sensitive", "Pyridine odor"]
  },
  "alternatives": ["TMO (Trimethyloxonium tetrafluoroborate)", "BF3-Methanol", "MSTFA"]
}
```

## Cherry Studio Config Key: `"--tools", "DerivatizationReagentSelector"`
