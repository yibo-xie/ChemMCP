# Tool #93: NuclearEquationBalance

## 基本信息
- **工具名称**: NuclearEquationBalance
- **版本**: 0.1.0
- **模块文件**: src/chemmcp/tools/nuclear_equation_balance.py
- **分类**: General
- **描述**: 核反应方程式配平（质量数+电荷数守恒）

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
{"code_input": {"reactants": "U-235 n", "products": "Ba-141 Kr-92 3n"}, "text_input": {"equation_str": "U-235 + n → Ba-141 + Kr-92 + 3n"}, "output": {"balanced": true, "total_mass_reactants": 236, "total_mass_products": 236, "total_charge_reactants": 92, "total_charge_products": 92}}
```

## Cherry Studio 配置

```json
{
  "mcpServers": {
    "ChemMCP_93_NuclearEquationBalance": {
      "command": "uv",
      "args": [
        "--directory", "/home/wave/ChemMCP",
        "run", "-m", "chemmcp",
        "--tools", "NuclearEquationBalance"
      ]
    }
  }
}
```

## 测试状态
✅ 所有测试通过
