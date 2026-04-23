# #240 TunnelingProbability

> **量子隧穿概率计算**
> Generated: 2026-04-23 22:06:08

## Tool Info

| Field | Value |
|-------|-------|
| Name | TunnelingProbability |
| Version | 0.1.0 |
| Description | Calculate quantum tunneling (barrier penetration) probability using WKB approximation and exact solutions for rectangular, triangular, Gaussian, Eckart barriers, and alpha decay. |
| Category | General |
| Tags | Quantum Mechanics, Tunneling, WKB Approximation, Barrier Penetration, Alpha Decay |

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
        "TunnelingProbability"
      ]
    }
  }
}
```

## Test Results

**All tests passed** (`test_quantum_chemistry_tools_231_240.py`) - 10/10 tools verified.

### Verified Outputs

- Rectangular (1nm): T=**1.5253e-02**, R=**0.9847**
- Alpha decay: T=**6.3065e-16**, G=**17.5**
- FN emission (1e10 V/m): T=**4.8288e-04**

## Usage Examples

### Code Interface
```python
from chemmcp.tools import TunnelingProbability
tool = TunnelingProbability()
result = tool.run_code(...)
```

### Text Interface
```
> TunnelingProbability <parameters>
```

## Notes
- Physical constants from CODATA recommended values
- All energies in eV unless otherwise noted
- Distances in Bohr radii (a0) or Angstroms (A)
