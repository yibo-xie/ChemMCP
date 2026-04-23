# #236 HuckelMethod

> **休克尔分子轨道法计算π电子体系**
> Generated: 2026-04-23 22:06:08

## Tool Info

| Field | Value |
|-------|-------|
| Name | HuckelMethod |
| Version | 0.1.0 |
| Description | Hückel MO theory for conjugated π-electron systems: linear/cyclic polyenes, aromatic compounds, and radicals. |
| Category | General |
| Tags | Quantum Chemistry, Hückel Method, Pi Electrons, Conjugated Systems, Aromaticity |

## MCP Configuration (Cherry Studio)

```json
{
  "mcpServers": {
    "ChemMCP": {
      "command": "uv",
      "args": [
        "--directory",
        "/home/wave/ChemMCP",
        "run",
        "-m",
        "chemmcp",
        "--tools",
        "HuckelMethod"
      ]
    }
  }
}
```

## Test Results

**All tests passed** (`test_quantum_chemistry_tools_231_240.py`) - 10/10 tools verified.

### Verified Outputs

- Ethene: E_pi = **2.000beta**, deloc = **0.000beta**
- Butadiene: E_pi = **4.472beta**, deloc = **0.472beta**
- Benzene: cyclic topology, orbital energies +/-2beta, +/-1beta verified

## Usage Examples

### Code Interface
```python
from chemmcp.tools import HuckelMethod
tool = HuckelMethod()
result = tool.run_code(...)
```

### Text Interface
```
> HuckelMethod <parameters>
```

## Notes
- Physical constants from CODATA recommended values
- All energies in eV unless otherwise noted
- Distances in Bohr radii (a0) or Angstroms (A)
