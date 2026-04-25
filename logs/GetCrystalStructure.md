# GetCrystalStructure

> **Version:** 0.1.0 | **Function:** `get_crystal_structure`

## Description
Query crystal structure type for elements or common ionic compounds.

## Implementation Details
Uses a built-in crystallographic database with structure type, coordination number, packing efficiency, lattice parameters, space group information for elements and common compounds.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('material` | str | N/A | Element symbol or compound formula (e.g., Fe, NaCl, CsCl, ZnS, CaF2)'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `material` | str | 'Material identifier' |
| `crystal_system` | str | 'Crystal system/structure type' |
| `packing_efficiency` | float | 'Atomic packing factor (if applicable |
| `lattice_parameters` | dict | 'Lattice constants' |
| `description` | str | 'Detailed structural description' |

## Examples
```python
examples = [
{'code_input': {'material': 'NaCl'}, 'text_input': {'material': 'NaCl'}, 'output': {'material': 'NaCl', 'crystal_system': 'rock salt (fcc)', 'coordination_number': [6, 6
]
```

## Cherry Studio JSON Config
```json
{
  "name": "GetCrystalStructure",
  "description": "Query crystal structure type for elements or common ionic compounds.",
  "function": "get_crystal_structure",
  "parameters": {
    "('material": {
      "type": "str",
      "description": "Element symbol or compound formula (e.g., Fe, NaCl, CsCl, ZnS, CaF2)'),"
    }
  }
}
```
