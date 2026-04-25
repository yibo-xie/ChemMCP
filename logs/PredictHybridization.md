# PredictHybridization

> **Version:** 0.1.0 | **Function:** `predict_hybridization`

## Description
Predict hybridization state of atoms in a molecule from its SMILES string.

## Implementation Details
Uses RDKit to analyze molecular structure and determines hybridization for each atom based on steric number (number of sigma bonds + lone pairs). Maps steric number to hybridization type (sp, sp2, sp3, sp3d, sp3d2).

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('smiles` | str | N/A | SMILES string of the molecule'), |
| `('atom_index` | int | N/A | Optional: specific atom index (0-based). If not provided, returns all atoms.'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `atom_hybridizations` | list | 'Hybridization for each atom with symbol and steric number' |
| `summary` | dict | 'Summary of unique hybridizations found' |

## Examples
```python
examples = [
{'code_input': {'smiles': 'CCO', 'atom_index': None}, 'text_input': {'smiles': 'ethanol'}, 'output': {'atom_hybridizations': [...
]
```

## Cherry Studio JSON Config
```json
{
  "name": "PredictHybridization",
  "description": "Predict hybridization state of atoms in a molecule from its SMILES string.",
  "function": "predict_hybridization",
  "parameters": {
    "('smiles": {
      "type": "str",
      "description": "SMILES string of the molecule'),"
    },
    "('atom_index": {
      "type": "int",
      "description": "Optional: specific atom index (0-based). If not provided, returns all atoms.'),"
    }
  }
}
```
