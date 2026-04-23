# #239 SelectionRulesChecker

> **光谱跃迁选择定则验证**
> Generated: 2026-04-23 22:06:08

## Tool Info

| Field | Value |
|-------|-------|
| Name | SelectionRulesChecker |
| Version | 0.1.0 |
| Description | Check spectroscopic transition selection rules for electric dipole, magnetic dipole, electric quadrupole, vibrational, rotational, Raman, electronic atomic/molecular transitions. |
| Category | General |
| Tags | Quantum Mechanics, Selection Rules, Spectroscopy, Transitions, Electric Dipole |

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
        "SelectionRulesChecker"
      ]
    }
  }
}
```

## Test Results

**All tests passed** (`test_quantum_chemistry_tools_231_240.py`) - 10/10 tools verified.

### Verified Outputs

- 1s->2p E1: **ALLOWED** (strong)
- 1s->2s E1: **FORBIDDEN** (Δl = +0 ✗ (requires Δl = ±1))
- Same-l M1: **ALLOWED**

## Usage Examples

### Code Interface
```python
from chemmcp.tools import SelectionRulesChecker
tool = SelectionRulesChecker()
result = tool.run_code(...)
```

### Text Interface
```
> SelectionRulesChecker <parameters>
```

## Notes
- Physical constants from CODATA recommended values
- All energies in eV unless otherwise noted
- Distances in Bohr radii (a0) or Angstroms (A)
