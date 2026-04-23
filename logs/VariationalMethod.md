# #233 VariationalMethod

> **变分法求解近似基态能量**
> Generated: 2026-04-23 22:06:08

## Tool Info

| Field | Value |
|-------|-------|
| Name | VariationalMethod |
| Version | 0.1.0 |
| Description | Approximate ground state energy using variational method with various trial wavefunctions for quantum systems. |
| Category | General |
| Tags | Quantum Mechanics, Variational Method, Approximation, Ground State, Optimization |

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
        "VariationalMethod"
      ]
    }
  }
}
```

## Test Results

**All tests passed** (`test_quantum_chemistry_tools_231_240.py`) - 10/10 tools verified.

### Verified Outputs

- Harmonic/Gaussian trial: E_var = **1.090437 eV**, E_exact = **1.090437 eV**
  - Relative error: **0.0000%**

## Usage Examples

### Code Interface
```python
from chemmcp.tools import VariationalMethod
tool = VariationalMethod()
result = tool.run_code(...)
```

### Text Interface
```
> VariationalMethod <parameters>
```

## Notes
- Physical constants from CODATA recommended values
- All energies in eV unless otherwise noted
- Distances in Bohr radii (a0) or Angstroms (A)
