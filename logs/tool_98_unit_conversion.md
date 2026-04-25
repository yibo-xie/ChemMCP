# Tool #98: UnitConversion

## 基本信息
- **工具名称**: UnitConversion
- **版本**: 0.1.0
- **模块文件**: src/chemmcp/tools/unit_conversion.py
- **分类**: General
- **描述**: 化学常用单位换算（温度、压强、能量、长度等）

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
{"code_input": {"value": 1.0, "from_unit": "atm", "to_unit": "Pa"}, "text_input": {"conversion_str": "1 atm to Pa"}, "output": {"value": 1.0, "from_unit": "atm", "to_unit": "Pa", "result": 101325.0, "result_unit": "Pa"}}
```

## Cherry Studio 配置

```json
{
  "mcpServers": {
    "ChemMCP_98_UnitConversion": {
      "command": "uv",
      "args": [
        "--directory", "/home/wave/ChemMCP",
        "run", "-m", "chemmcp",
        "--tools", "UnitConversion"
      ]
    }
  }
}
```

## 测试状态
✅ 所有测试通过
