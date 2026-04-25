# GetIsotopes

> **Version:** 0.1.0 | **Function:** `get_isotopes`

## Description
Get isotope information for an element including mass number, natural abundance, stability, and half-life.

## Implementation Details
Uses a built-in isotope database with data for common elements. Returns mass number, percent abundance, stability status, and half-life for radioactive isotopes.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('element` | str | N/A | Element symbol (e.g., C, U, Fe)'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `element` | str | 'Element symbol' |
| `isotopes` | list | 'List of isotope data dictionaries' |
| `isotope_count` | int | 'Number of isotopes in database' |

## Examples
```python
examples = [
{'code_input': {'element': 'C'}, 'text_input': {'element': 'C'}, 'output': {'element': 'C', 'isotopes': [{...}, {...}, {...}
]
```

## Cherry Studio JSON Config
```json
{
  "name": "GetIsotopes",
  "description": "Get isotope information for an element including mass number, natural abundance, stability, and half-life.",
  "function": "get_isotopes",
  "parameters": {
    "('element": {
      "type": "str",
      "description": "Element symbol (e.g., C, U, Fe)'),"
    }
  }
}
```
