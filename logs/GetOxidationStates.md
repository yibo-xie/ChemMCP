# GetOxidationStates

> **Version:** 0.1.0 | **Function:** `get_oxidation_states`

## Description
Query common oxidation states of an element with stability information.

## Implementation Details
Uses a built-in database of oxidation states for all common elements. Returns all known oxidation states, the most common ones, and descriptive notes for each state including typical compounds.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('element` | str | N/A | Element symbol (e.g., Fe, Mn, Cl)'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `element` | str | 'Element symbol' |
| `oxidation_states` | list | 'All known oxidation states' |
| `state_details` | dict | 'Details for each oxidation state with examples' |

## Examples
```python
examples = [
{'code_input': {'element': 'Fe'}, 'text_input': {'element': 'Fe'}, 'output': {'element': 'Fe', 'oxidation_states': [+2, +3, +6
]
```

## Cherry Studio JSON Config
```json
{
  "name": "GetOxidationStates",
  "description": "Query common oxidation states of an element with stability information.",
  "function": "get_oxidation_states",
  "parameters": {
    "('element": {
      "type": "str",
      "description": "Element symbol (e.g., Fe, Mn, Cl)'),"
    }
  }
}
```
