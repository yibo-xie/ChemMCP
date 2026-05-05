# Log #303: ExtractionOptimizer (液-液萃取优化)

## Tool Info
- **Tool ID**: 303
- **Class Name**: `ExtractionOptimizer`
- **Module**: `extraction_optimizer`
- **Version**: 0.1.0

## Core Logic
液-液萃取条件优化，基于分配系数 (Kd/Kow) 计算单次和多次萃取效率。支持 pH 依赖性优化（根据 pKa 计算中性分子比例），自动推荐最佳萃取次数和溶剂体积比。

## Input Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| analyte | str | N/A | 分析物名称/标识 |
| pKa | float | N/A | 酸解离常数 |
| logP | float | N/A | 油水分配系数 |
| extraction_solvent | str | "diethyl_ether" | 萃取溶剂 |
| aqueous_phase_pH | float | 7.0 | 水相 pH |
| num_extractions | int | 1 | 萃取次数 |
| solvent_volume_ratio | float | 1.0 | 有机相:水相体积比 |

## Output Fields
| Field | Type | Description |
|-------|------|-------------|
| analyte | str | 分析物名称 |
| distribution_coefficient_Kd | float | 分配系数 Kd |
| extraction_efficiency_pct | float | 单次萃取效率 (%) |
| total_efficiency_pct | float | 总萃取效率 (%) |
| ph_optimization | dict | pH 优化建议 |
| recommendations | list | 萃取条件优化建议 |

## Example: Benzoic Acid / Diethyl Ether / 1x
```json
{
  "analyte": "benzoic_acid",
  "distribution_coefficient_Kd": 74.13,
  "extraction_efficiency_pct": 97.4,
  "total_efficiency_pct": 97.4,
  "ph_optimization": {
    "optimal_pH_range": "pH < 2.0 (2 units below pKa)",
    "neutral_fraction_at_ph7": 0.0016,
    "recommendation": "Acidify to pH 2-3 for efficient extraction"
  }
}
```

## Cherry Studio Config Key: `"--tools", "ExtractionOptimizer"`
