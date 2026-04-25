# GetElectronConfiguration

> **Version:** 0.1.0 | **Function:** `get_electron_configuration`

## Description
Get the electron configuration of an element, including full and noble-gas shorthand forms.

## Implementation Details
Uses the built-in periodic table database to return both the full electron configuration (e.g., '1s² 2s² 2p⁶ 3s² 3p⁴') and the noble-gas shorthand (e.g., '[Ne] 3s² 3p⁴').

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('element` | str | N/A | Element symbol (e.g., S, Fe)'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `full_config` | str | 'Full electron configuration' |
| `noble_config` | str | 'Noble gas shorthand notation' |
| `symbol` | str | 'Element symbol' |

## Examples
```python
examples = [
{'code_input': {'element': 'S'}, 'text_input': {'element': 'S'}, 'output': {'full_config': '[Ne
]
```

## Cherry Studio JSON Config
```json
{
  "name": "GetElectronConfiguration",
  "description": "Get the electron configuration of an element, including full and noble-gas shorthand forms.",
  "function": "get_electron_configuration",
  "parameters": {
    "('element": {
      "type": "str",
      "description": "Element symbol (e.g., S, Fe)'),"
    }
  }
}
```
