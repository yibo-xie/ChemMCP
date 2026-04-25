# GetElectronAffinity

> **Version:** 0.1.0 | **Function:** `get_electron_affinity`

## Description
Query electron affinity of an element in kJ/mol.

## Implementation Details
Returns the first electron affinity value from NIST/CRC Handbook data. Negative values indicate energy release (exothermic). Returns None for elements that do not readily accept an electron.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('element` | str | N/A | Element symbol (e.g., Cl, O)'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `element` | str | 'Element symbol' |
| `electron_affinity_kj_mol` | float | 'Electron affinity in kJ/mol' |
| `note` | str | 'Interpretation note' |

## Examples
```python
examples = [
{'code_input': {'element': 'Cl'}, 'text_input': {'element': 'Cl'}, 'output': {'element': 'Cl', 'electron_affinity_kj_mol': 349.0, 'note': '...'}},
        {'code_input': {'element': 'N'}, 'text_input': {'element': 'N'}, 'output': {'element': 'N', 'electron_affinity_kj_mol': -6.8, 'note': '...'}},
]
```

## Cherry Studio JSON Config
```json
{
  "name": "GetElectronAffinity",
  "description": "Query electron affinity of an element in kJ/mol.",
  "function": "get_electron_affinity",
  "parameters": {
    "('element": {
      "type": "str",
      "description": "Element symbol (e.g., Cl, O)'),"
    }
  }
}
```
