# #235 MolecularOrbitalDiagram

> **分子轨道能级图生成**
> Generated: 2026-04-23 22:06:08

## Tool Info

| Field | Value |
|-------|-------|
| Name | MolecularOrbitalDiagram |
| Version | 0.1.0 |
| Description | Generate molecular orbital diagrams for diatomic and simple polyatomic molecules using LCAO method. |
| Category | General |
| Tags | Quantum Chemistry, Molecular Orbitals, LCAO, Bond Order, HOMO-LUMO |

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
        "MolecularOrbitalDiagram"
      ]
    }
  }
}
```

## Test Results

**All tests passed** (`test_quantum_chemistry_tools_231_240.py`) - 10/10 tools verified.

### Verified Outputs

- O2: BO=**2.0**, **paramagnetic**, gap=**8.5 eV**
- N2: BO=**3.0**, **diamagnetic**, gap=**22.0 eV**
- H2O: **bent, 104.5°**, **C2v**, gap=**18.5 eV**

## Usage Examples

### Code Interface
```python
from chemmcp.tools import MolecularOrbitalDiagram
tool = MolecularOrbitalDiagram()
result = tool.run_code(...)
```

### Text Interface
```
> MolecularOrbitalDiagram <parameters>
```

## Notes
- Physical constants from CODATA recommended values
- All energies in eV unless otherwise noted
- Distances in Bohr radii (a0) or Angstroms (A)
