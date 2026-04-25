# Tool #92: HalfLifeCalculation

## 基本信息
- **工具名称**: HalfLifeCalculation
- **版本**: 0.1.0
- **模块文件**: src/chemmcp/tools/half_life_calculation.py
- **分类**: General
- **描述**: 半衰期相关计算（剩余量、衰变常数、经历时间、初始量）

## 输入输出签名

### Code Input (_run_base)
| 参数名 | 类型 | 默认值 | 描述 |
|--------|------|--------|------|
(见源码 code_input_sig)

### Text Input (_run_text)
| 参数名 | 类型 | 默认值 | 描述 |
|--------|------|--------|------|
(见源码 text_input_sig)

### Output
| 字段名 | 类型 | 描述 |
|--------|------|------|
(见源码 output_sig)

## 使用示例

```json
{"code_input": {"calc_type": "remaining_amount", "half_life": 5730.0, "initial_amount": 100.0, "time_elapsed": 17190.0}, "text_input": {"query_str": "C-14 after 3 half-lives"}, "output": {"calc_type": "remaining_amount", "half_life": 5730.0, "initial_amount": 100.0, "time_elapsed": 17190.0, "decay_constant": 0.0001209681, "remaining_amount": 12.5, "half_lives_elapsed": 3.0}}
```

## Cherry Studio 配置

```json
{
  "mcpServers": {
    "ChemMCP_92_HalfLifeCalculation": {
      "command": "uv",
      "args": [
        "--directory", "/home/wave/ChemMCP",
        "run", "-m", "chemmcp",
        "--tools", "HalfLifeCalculation"
      ]
    }
  }
}
```

## 测试状态
✅ 所有测试通过
