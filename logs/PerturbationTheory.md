# #234 PerturbationTheory

> **微扰理论能量修正计算**
> Generated: 2026-04-23 22:06:08

## Tool Info

| Field | Value |
|-------|-------|
| Name | PerturbationTheory |
| Version | 0.1.0 |
| Description | Calculate energy corrections using non-degenerate and degenerate perturbation theory for quantum systems. |
| Category | General |
| Tags | Quantum Mechanics, Perturbation Theory, Energy Correction, Approximation, Stark Effect |

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
        "PerturbationTheory"
      ]
    }
  }
}
```

## Test Results

**All tests passed** (`test_quantum_chemistry_tools_231_240.py`) - 10/10 tools verified.

### Verified Outputs

- Two-level system: E_PT = **-3.465366e-01 eV**, error = **3.3152%**
- He ground state: E_PT = **-74.83 eV**, error ~**5.3%**

## Usage Examples

### Code Interface
```python
from chemmcp.tools import PerturbationTheory
tool = PerturbationTheory()
result = tool.run_code(...)
```

### Text Interface
```
> PerturbationTheory <parameters>
```

## Notes
- Physical constants from CODATA recommended values
- All energies in eV unless otherwise noted
- Distances in Bohr radii (a0) or Angstroms (A)
