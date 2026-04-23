# #238 SpinOrbitCoupling

> **自旋-轨道耦合能级计算**
> Generated: 2026-04-23 22:06:08

## Tool Info

| Field | Value |
|-------|-------|
| Name | SpinOrbitCoupling |
| Version | 0.1.0 |
| Description | Calculate spin-orbit coupling energy level splitting for atoms and fine structure analysis. |
| Category | General |
| Tags | Quantum Mechanics, Spin-Orbit Coupling, Fine Structure, Atomic Physics, Term Symbols |

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
        "SpinOrbitCoupling"
      ]
    }
  }
}
```

## Test Results

**All tests passed** (`test_quantum_chemistry_tools_231_240.py`) - 10/10 tools verified.

### Verified Outputs

- Na 3p: zeta=**17.20 cm-1**, dE=**25.79 cm-1**
  - Terms: **['2P_{0.5}', '2P_{1.5}']**
- Na 3s: No splitting (**1 level**) OK

## Usage Examples

### Code Interface
```python
from chemmcp.tools import SpinOrbitCoupling
tool = SpinOrbitCoupling()
result = tool.run_code(...)
```

### Text Interface
```
> SpinOrbitCoupling <parameters>
```

## Notes
- Physical constants from CODATA recommended values
- All energies in eV unless otherwise noted
- Distances in Bohr radii (a0) or Angstroms (A)
