# Tool #94: BindingEnergy

## 基本信息
- **工具名称**: BindingEnergy
- **版本**: 0.1.0
- **模块文件**: src/chemmcp/tools/binding_energy.py
- **分类**: General
- **描述**: 计算核结合能及比结合能（MeV/nucleon）

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
{"code_input": {"nuclide": "Fe-56", "custom_mass_amu": 0.0, "output_unit": "MeV"}, "text_input": {"nuclide_str": "Fe-56"}, "output": {"nuclide": "Fe-56", "mass_number": 56, "atomic_number": 26, "binding_energy_per_nucleon_mev": 8.79}}
```

## Cherry Studio 配置

```json
{
  "mcpServers": {
    "ChemMCP_94_BindingEnergy": {
      "command": "uv",
      "args": [
        "--directory", "/home/wave/ChemMCP",
        "run", "-m", "chemmcp",
        "--tools", "BindingEnergy"
      ]
    }
  }
}
```

## 测试状态
✅ 所有测试通过
