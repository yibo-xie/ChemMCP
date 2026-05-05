# Log #302: StandardSolutionPrep (标准溶液配制)

## Tool Info
- **Tool ID**: 302
- **Class Name**: `StandardSolutionPrep`
- **Module**: `standard_solution_prep`
- **Version**: 0.1.0

## Core Logic
标准溶液配制指导，包括称量计算、定容步骤。支持摩尔浓度 (mol/L) 和质量浓度 (g/L) 两种模式。自动计算所需溶质质量，考虑纯度校正因子，生成完整的分步配制流程。

## Input Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| solute | str | N/A | 溶质名称或化学式（如 K2Cr2O7, NaCl） |
| concentration | float | N/A | 目标浓度值 |
| concentration_unit | str | "mol/L" | 浓度单位 (mol/L 或 g/L) |
| final_volume_ml | float | N/A | 定容体积 (mL) |
| purity_pct | float | 100.0 | 溶质纯度 (%) |

## Output Fields
| Field | Type | Description |
|-------|------|-------------|
| solute_info | dict | 溶质信息（名称、分子式、MW） |
| mass_required_g | float | 需称取的质量 (g)，已考虑纯度 |
| purity_correction | float | 纯度校正因子 |
| preparation_steps | list | 分步操作指南 |
| notes | list | 注意事项和提示 |

## Example: K2Cr2O7 0.02M / 250mL
```json
{
  "solute_info": {"name": "Potassium dichromate", "formula": "K2Cr2O7", "molar_mass": 294.18},
  "mass_required_g": 1.4709,
  "purity_correction": 1.0,
  "preparation_steps": [
    {"step": 1, "action": "Calculate required mass", "detail": "1.4709 g of K2Cr2O7"},
    {"step": 2, "action": "Weigh solute", "detail": "Using analytical balance (±0.0001 g)"},
    ...
  ],
  "notes": ["K2Cr2O7 is a strong oxidizer - handle with gloves", "Store in amber bottle"]
}
```

## Cherry Studio Config Key: `"--tools", "StandardSolutionPrep"`
