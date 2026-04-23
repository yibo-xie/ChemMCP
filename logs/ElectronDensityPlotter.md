# #237 ElectronDensityPlotter

> **电子密度分布可视化**
> Generated: 2026-04-23 22:06:08

## Tool Info

| Field | Value |
|-------|-------|
| Name | ElectronDensityPlotter |
| Version | 0.1.0 |
| Description | Calculate and visualize electron density distribution for hydrogen-like atoms and quantum systems. |
| Category | General |
| Tags | Quantum Mechanics, Electron Density, Visualization, Atomic Orbitals, Radial Distribution |

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
        "ElectronDensityPlotter"
      ]
    }
  }
}
```

## Test Results

**All tests passed** (`test_quantum_chemistry_tools_231_240.py`) - 10/10 tools verified.

### Verified Outputs

- H 1s: r_mp(RDF)=**1.003 a0**, <r>=**1.497 a0**
- H 2p: r_mp(RDF)=**4.022 a0**, radial_nodes=**0**

## Usage Examples

### Code Interface
```python
from chemmcp.tools import ElectronDensityPlotter
tool = ElectronDensityPlotter()
result = tool.run_code(...)
```

### Text Interface
```
> ElectronDensityPlotter <parameters>
```

## Notes
- Physical constants from CODATA recommended values
- All energies in eV unless otherwise noted
- Distances in Bohr radii (a0) or Angstroms (A)
