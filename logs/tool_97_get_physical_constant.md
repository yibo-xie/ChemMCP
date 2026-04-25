# Tool #97: GetPhysicalConstant

## 基本信息
- **工具名称**: GetPhysicalConstant
- **版本**: 0.1.0
- **模块文件**: src/chemmcp/tools/get_physical_constant.py
- **分类**: General
- **描述**: 查询物理常数（NA、R、F、h、c等25个常数）

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
{"code_input": {"constant_name": "Avogadro"}, "text_input": {"query_str": "Avogadro number"}, "output": {"constant_name": "Avogadro constant", "symbol": "Nₐ", "value": 6.02214076e+23, "unit": "mol⁻¹"}}
```

## Cherry Studio 配置

```json
{
  "mcpServers": {
    "ChemMCP_97_GetPhysicalConstant": {
      "command": "uv",
      "args": [
        "--directory", "/home/wave/ChemMCP",
        "run", "-m", "chemmcp",
        "--tools", "GetPhysicalConstant"
      ]
    }
  }
}
```

## 测试状态
✅ 所有测试通过
