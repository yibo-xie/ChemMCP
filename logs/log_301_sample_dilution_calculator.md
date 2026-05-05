# Log #301: SampleDilutionCalculator (样品稀释计算器)

## Tool Info
- **Tool ID**: 301
- **Class Name**: `SampleDilutionCalculator`
- **Module**: `sample_dilution_calculator`
- **Version**: 0.1.0
- **Category**: Sample Preparation & Analytical Chemistry

## Core Logic
计算稀释比例、终浓度和所需溶剂体积。支持单步稀释和多步 serial dilution（连续稀释）。基于 C₁V₁ = C₂V₂ 公式，自动计算稀释因子、需添加的溶剂体积、终浓度，并生成每步操作的详细步骤。

## Input Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| initial_concentration | float | N/A | 初始浓度（与目标浓度单位一致） |
| initial_volume | float | N/A | 初始体积 (mL) |
| final_volume | float | N/A | 稀释后总体积 (mL) |
| dilution_steps | int | 1 | 连续稀释步数（1=单步稀释） |

## Output Fields
| Field | Type | Description |
|-------|------|-------------|
| dilution_factor | float | 总稀释倍数 |
| solvent_volume_needed | float | 需添加的溶剂体积 (mL) |
| final_concentration | float | 稀释后终浓度 |
| dilution_ratio | str | 稀释比表示（如 "1:100"） |
| step_details | list | 每步操作详情（step, dilution_factor, solvent_added_ml, conc） |

## Example Usage

### Code Input
```python
tool = SampleDilutionCalculator()
result = tool.run_code(
    initial_concentration=1000.0,
    initial_volume=1.0,
    final_volume=100.0,
    dilution_steps=1
)
```

### Output
```json
{
  "dilution_factor": 100.0,
  "solvent_volume_needed": 99.0,
  "final_concentration": 10.0,
  "dilution_ratio": "1:100",
  "step_details": [
    {"step": 1, "dilution_factor": 100.0, "solvent_added_ml": 99.0, "conc": 10.0}
  ]
}
```

### Text Interface
```python
result = tool.run_text("500 2.0 50.0")
# → dilution_factor=25.0, final_conc=20.0
```

## Key Implementation Notes
- 支持 serial dilution：自动将总稀释因子分解为 N 步等比稀释
- 内置浓度/体积验证：终浓度不能超过初始浓度
- 使用 `logging.getLogger(__name__)` 记录计算过程
- 无需外部 API key（纯计算工具）
- 错误处理：无效输入抛出 `ChemMCPError`

## Cherry Studio JSON Config
```json
{
  "mcpServers": {
    "ChemMCP_301_SampleDilutionCalculator": {
      "command": "uv",
      "args": ["--directory", "/home/wave/ChemMCP", "run", "-m", "chemmcp", "--tools", "SampleDilutionCalculator"],
      "env": {}
    }
  }
}
```
