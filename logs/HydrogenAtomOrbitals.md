# #231 HydrogenAtomOrbitals

> **氢原子轨道可视化和能级计算**
> Generated: 2026-04-23 22:06:08

## Tool Info

| Field | Value |
|-------|-------|
| Name | HydrogenAtomOrbitals |
| Version | 0.1.0 |
| Description | Calculate hydrogen atom orbital properties: energy levels, radial wavefunction, probability density, orbital shapes, and node structure. |
| Category | General |
| Tags | Quantum Mechanics, Hydrogen Atom, Orbitals, Energy Levels, Wavefunction |

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
        "HydrogenAtomOrbitals"
      ]
    }
  }
}
```

## Test Results

**All tests passed** (`test_quantum_chemistry_tools_231_240.py`) - 10/10 tools verified.

### Verified Outputs

- Ground state (1s): E = **-13.6057 eV**, r_mp = **1.000 a0**, <r> = **1.500 a0**
- Excited state (2p): E = **-3.4014 eV**, nodes = **1**
- He+ ion (Z=2): E = **-54.42 eV**, r_mp = **0.500 a0**

## Usage Examples

### Code Interface
```python
from chemmcp.tools import HydrogenAtomOrbitals
tool = HydrogenAtomOrbitals()
result = tool.run_code(...)
```

### Text Interface
```
> HydrogenAtomOrbitals <parameters>
```

## Notes
- Physical constants from CODATA recommended values
- All energies in eV unless otherwise noted
- Distances in Bohr radii (a0) or Angstroms (A)
