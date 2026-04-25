# GetIonizationEnergy

> **Version:** 0.1.0 | **Function:** `get_ionization_energy`

## Description
Get ionization energy data for an element (first through available ionizations).

## Implementation Details
Uses NIST-standard ionization energy data in kJ/mol. Returns all available successive ionization energies from IE1 upward.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('element` | str | N/A | Element symbol (e.g., Na, Fe)'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `element` | str | 'Element symbol' |
| `ionization_energies` | dict | 'IE1, IE2, ... in kJ/mol' |
| `unit` | str | 'Unit of measurement' |

## Examples
```python
examples = [
{'code_input': {'element': 'Na'}, 'text_input': {'element': 'Na'}, 'output': {'element': 'Na', 'ionization_energies': {'IE1': 496, 'IE2': 4563}, 'unit': 'kJ/mol'}},
]
```

## Cherry Studio JSON Config
```json
{
  "name": "GetIonizationEnergy",
  "description": "Get ionization energy data for an element (first through available ionizations).",
  "function": "get_ionization_energy",
  "parameters": {
    "('element": {
      "type": "str",
      "description": "Element symbol (e.g., Na, Fe)'),"
    }
  }
}
```
