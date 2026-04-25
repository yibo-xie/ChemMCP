# Tool #99: SignificantFigures

## 基本信息
- **工具名称**: SignificantFigures
- **版本**: 0.1.0
- **模块文件**: src/chemmcp/tools/significant_figures.py
- **分类**: General
- **描述**: 有效数字处理（计数、四舍五入、科学计数法、运算）

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
{"code_input": {"operation": "count", "value_a": "0.004500"}, "text_input": {"query_str": "count sig figs of 0.004500"}, "output": {"operation": "count", "value_a": "0.004500", "sig_figs": 4}}
```

## Cherry Studio 配置

```json
{
  "mcpServers": {
    "ChemMCP_99_SignificantFigures": {
      "command": "uv",
      "args": [
        "--directory", "/home/wave/ChemMCP",
        "run", "-m", "chemmcp",
        "--tools", "SignificantFigures"
      ]
    }
  }
}
```

## 测试状态
✅ 所有测试通过
