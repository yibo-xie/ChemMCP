# #232 SchrodingerSolver1d

> **一维薛定谔方程数值求解**
> Generated: 2026-04-23 22:06:08

## Tool Info

| Field | Value |
|-------|-------|
| Name | SchrodingerSolver1d |
| Version | 0.1.0 |
| Description | Numerically solve 1D time-independent Schrödinger equation for various potential wells using finite difference method. |
| Category | General |
| Tags | Quantum Mechanics, Schrodinger Equation, Numerical Methods, Eigenvalues, Finite Difference |

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
        "SchrodingerSolver1d"
      ]
    }
  }
}
```

## Test Results

**All tests passed** (`test_quantum_chemistry_tools_231_240.py`) - 10/10 tools verified.

### Verified Outputs

- Infinite square well: E0 = **0.3760 eV**, E1 = **1.5041 eV**
  - Energy ratio E1/E0 = **4.00** (expected ~4.00)

## Usage Examples

### Code Interface
```python
from chemmcp.tools import SchrodingerSolver1d
tool = SchrodingerSolver1d()
result = tool.run_code(...)
```

### Text Interface
```
> SchrodingerSolver1d <parameters>
```

## Notes
- Physical constants from CODATA recommended values
- All energies in eV unless otherwise noted
- Distances in Bohr radii (a0) or Angstroms (A)
