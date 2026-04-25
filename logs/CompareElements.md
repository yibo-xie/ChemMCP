# CompareElements

> **Version:** 0.1.0 | **Function:** `compare_elements`

## Description
Compare properties of multiple elements side by side.

## Implementation Details
Accepts a list of element symbols and a property name, returns a comparison table with values sorted by rank. Supports atomic weight, electronegativity, atomic number, first ionization energy, electron affinity, period, and group.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `comparison` | dict | 'Ranked comparison table with element values and ranking' |
| `property` | str | 'Property that was compared' |
| `trend_note` | str | 'Brief trend analysis' |

## Examples
```python
examples = [
{'code_input': {'elements': ['Li', 'Na', 'K'
]
```

## Cherry Studio JSON Config
```json
{
  "name": "CompareElements",
  "description": "Compare properties of multiple elements side by side.",
  "function": "compare_elements",
  "parameters": {}
}
```
