# GetElementInfo

> **Version:** 0.1.0 | **Function:** `get_element_info`

## Description
Get complete element information by element symbol or atomic number.

## Implementation Details
Uses a built-in periodic table database with IUPAC data for all 118 elements. Returns atomic number, symbol, name, atomic weight, Pauling electronegativity, electron configuration, category, group, period, and block.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('element` | str | N/A | Element symbol (e.g., Fe) or atomic number (e.g., 26)'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `info` | dict | 'Complete element information dictionary' |

## Examples
```python
examples = [
{'code_input': {'element': 'O'}, 'text_input': {'element': 'O'}, 'output': {'info': {'atomic_number': 8, 'symbol': 'O'}}},
        {'code_input': {'element': 26}, 'text_input': {'element': '26'}, 'output': {'info': {'atomic_number': 26, 'symbol': 'Fe'}}},
]
```

## Cherry Studio JSON Config
```json
{
  "name": "GetElementInfo",
  "description": "Get complete element information by element symbol or atomic number.",
  "function": "get_element_info",
  "parameters": {
    "('element": {
      "type": "str",
      "description": "Element symbol (e.g., Fe) or atomic number (e.g., 26)'),"
    }
  }
}
```
