# GetBondLength

> **Version:** 0.1.0 | **Function:** `get_bond_length`

## Description
Query standard bond lengths in picometers (pm) or Angstroms.

## Implementation Details
Uses a database of average experimental bond lengths from CRC Handbook and crystallographic data. Returns bond length for a given element pair and bond type (single, double, triple, aromatic). Can also estimate from covalent radii sum if exact data is unavailable.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('element1` | str | N/A | First element symbol (e.g., C)'), |
| `('element2` | str | N/A | Second element symbol (e.g., O)'), |
| `('bond_type` | str | single | Bond type: single, double, triple, or aromatic'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `bond` | str | 'Bond specification' |
| `length_pm` | float | 'Bond length in picometers' |
| `length_angstrom` | float | 'Bond length in Angstroms' |
| `source` | str | 'Data source (experimental/covalent radii estimate |

## Examples
```python
examples = [
{'code_input': {'element1': 'C', 'element2': 'C', 'bond_type': 'single'}, 'text_input': {'query': 'C-C single'}, 'output': {'length_pm': 154, 'length_angstrom': 1.54, 'bond': 'C-C single', 'source': 'CRC Handbook'}},
        {'code_input': {'element1': 'C', 'element2': 'O', 'bond_type': 'double'}, 'text_input': {'query': 'C=O double'}, 'output': {'length_pm': 123, 'length_angstrom': 1.23, 'bond': 'C=O double', 'source': 'CRC Handbook'}},
        {'code_input': {'element1': 'C', 'element2': 'N', 'bond_type': 'triple'}, 'text_input': {'query': 'C≡N triple'}, 'output': {'length_pm': 116, 'length_angstrom': 1.16, 'bond': 'C≡N triple', 'source': 'CRC Handbook'}},
]
```

## Cherry Studio JSON Config
```json
{
  "name": "GetBondLength",
  "description": "Query standard bond lengths in picometers (pm) or Angstroms.",
  "function": "get_bond_length",
  "parameters": {
    "('element1": {
      "type": "str",
      "description": "First element symbol (e.g., C)'),"
    },
    "('element2": {
      "type": "str",
      "description": "Second element symbol (e.g., O)'),"
    },
    "('bond_type": {
      "type": "str",
      "description": "Bond type: single, double, triple, or aromatic'),"
    }
  }
}
```
