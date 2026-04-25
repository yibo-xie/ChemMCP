# Tool #91: RadioactiveDecay

## 基本信息
- **工具名称**: RadioactiveDecay
- **版本**: 0.1.0
- **模块文件**: src/chemmcp/tools/radioactive_decay.py
- **分类**: General
- **描述**: 放射性衰变类型判断（α/β⁻/β⁺/γ/IT）

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
{"code_input": {"nuclide": "U-238"}, "text_input": {"nuclide_str": "U-238"}, "output": {"nuclide": "U-238", "decay_mode": "α", "daughter_nuclide": "Th-234", "half_life_years": 4468000000.0}}
```

## Cherry Studio 配置

```json
{
  "mcpServers": {
    "ChemMCP_91_RadioactiveDecay": {
      "command": "uv",
      "args": [
        "--directory", "/home/wave/ChemMCP",
        "run", "-m", "chemmcp",
        "--tools", "RadioactiveDecay"
      ]
    }
  }
}
```

## 测试状态
✅ 所有测试通过
