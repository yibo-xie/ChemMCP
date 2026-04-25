# Tool #95: MassDefect

## 基本信息
- **工具名称**: MassDefect
- **版本**: 0.1.0
- **模块文件**: src/chemmcp/tools/mass_defect.py
- **分类**: General
- **描述**: 计算核质量亏损及其能量当量

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
{"code_input": {"nuclide": "He-4", "custom_mass_amu": 0.0}, "text_input": {"nuclide_str": "He-4"}, "output": {"nuclide": "He-4", "proton_count": 2, "neutron_count": 2, "mass_defect_u": 0.030377, "energy_equivalent_mev": 28.30}}
```

## Cherry Studio 配置

```json
{
  "mcpServers": {
    "ChemMCP_95_MassDefect": {
      "command": "uv",
      "args": [
        "--directory", "/home/wave/ChemMCP",
        "run", "-m", "chemmcp",
        "--tools", "MassDefect"
      ]
    }
  }
}
```

## 测试状态
✅ 所有测试通过
