# Tool #100: DimensionalAnalysis

## 基本信息
- **工具名称**: DimensionalAnalysis
- **版本**: 0.1.0
- **模块文件**: src/chemmcp/tools/dimensional_analysis.py
- **分类**: General
- **描述**: 量纲分析辅助（查询量纲、检查一致性、推导量纲）

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
{"code_input": {"operation": "query", "quantity": "pressure"}, "text_input": {"query_str": "dimensions of pressure"}, "output": {"quantity": "pressure", "dimensional_formula": "M.L⁻¹.T⁻²"}}
```

## Cherry Studio 配置

```json
{
  "mcpServers": {
    "ChemMCP_100_DimensionalAnalysis": {
      "command": "uv",
      "args": [
        "--directory", "/home/wave/ChemMCP",
        "run", "-m", "chemmcp",
        "--tools", "DimensionalAnalysis"
      ]
    }
  }
}
```

## 测试状态
✅ 所有测试通过
