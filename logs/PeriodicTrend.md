# PeriodicTrend

> **Version:** 0.1.0 | **Function:** `periodic_trend`

## Description
Query periodic table trends for properties like atomic radius, electronegativity, ionization energy, electron affinity, and metallic character.

## Implementation Details
Returns trend descriptions (how a property changes across periods and down groups), plus actual data values for specific periods or groups when requested.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('property` | str | N/A | Property to query: atomic_radius, electronegativity, ionization_energy, electron_affinity, metallic_character'), |
| `('period` | int | N/A | Optional: show data for this period number (1-7)'), |
| `('group` | int | N/A | Optional: show data for this group number (1-18)'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `property` | str | 'Property queried' |
| `trend_description` | str | 'How the property changes across periods and groups' |
| `data` | dict | 'Actual values if period or group specified' |

## Examples
```python
examples = [
{'code_input': {'property': 'electronegativity', 'period': None, 'group': None}, 'text_input': {'query': 'electronegativity'}, 'output': {'property': 'electronegativity', 'trend_description': '...', 'data': {}}},
        {'code_input': {'property': 'atomic_radius', 'period': 3, 'group': None}, 'text_input': {'query': 'atomic radius period 3'}, 'output': {'property': 'atomic_radius', 'trend_description': '...', 'data': {}}},
]
```

## Cherry Studio JSON Config
```json
{
  "name": "PeriodicTrend",
  "description": "Query periodic table trends for properties like atomic radius, electronegativity, ionization energy, electron affinity, and metallic character.",
  "function": "periodic_trend",
  "parameters": {
    "('property": {
      "type": "str",
      "description": "Property to query: atomic_radius, electronegativity, ionization_energy, electron_affinity, metallic_character'),"
    },
    "('period": {
      "type": "int",
      "description": "Optional: show data for this period number (1-7)'),"
    },
    "('group": {
      "type": "int",
      "description": "Optional: show data for this group number (1-18)'),"
    }
  }
}
```
