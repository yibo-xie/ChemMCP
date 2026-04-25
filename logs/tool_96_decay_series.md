# Tool #96: DecaySeries

## 基本信息
- **工具名称**: DecaySeries
- **版本**: 0.1.0
- **模块文件**: src/chemmcp/tools/decay_series.py
- **分类**: General
- **描述**: 查询三大天然放射系（铀系、锕系、钍系）

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
{"code_input": {"series_name": "uranium"}, "text_input": {"query_str": "U-238 decay series"}, "output": {"series_name": "Uranium (4n+2)", "parent": "U-238", "final_stable": "Pb-206", "total_steps": 15}}
```

## Cherry Studio 配置

```json
{
  "mcpServers": {
    "ChemMCP_96_DecaySeries": {
      "command": "uv",
      "args": [
        "--directory", "/home/wave/ChemMCP",
        "run", "-m", "chemmcp",
        "--tools", "DecaySeries"
      ]
    }
  }
}
```

## 测试状态
✅ 所有测试通过
