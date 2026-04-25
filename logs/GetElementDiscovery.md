# GetElementDiscovery

> **Version:** 0.1.0 | **Function:** `get_element_discovery`

## Description
Get discovery history of an element including discoverer, year, place of discovery, name origin, and etymology.

## Implementation Details
Uses a built-in historical database covering elements discovered from antiquity through the modern era. Returns discoverer(s), year, location, naming origin, and interesting historical notes.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('element` | str | N/A | Element symbol (e.g., O, U, He)'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `element` | str | 'Element symbol' |
| `discoverer` | str | 'Discoverer(s |
| `place` | str | 'Place of discovery' |
| `name_origin` | str | 'Origin of the element name' |
| `etymology` | str | 'Etymological details' |

## Examples
```python
examples = [
{'code_input': {'element': 'O'}, 'text_input': {'element': 'O'}, 'output': {'element': 'O', 'discoverer': 'Carl Wilhelm Scheele', 'year': 1774, 'place': 'Sweden', 'name_origin': 'Greek', 'etymology': '...'}},
        {'code_input': {'element': 'U'}, 'text_input': {'element': 'U'}, 'output': {'element': 'U', 'discoverer': 'Martin Heinrich Klaproth', 'year': 1789, 'place': 'Germany', 'name_origin': 'Uranus', 'etymology': '...'}},
]
```

## Cherry Studio JSON Config
```json
{
  "name": "GetElementDiscovery",
  "description": "Get discovery history of an element including discoverer, year, place of discovery, name origin, and etymology.",
  "function": "get_element_discovery",
  "parameters": {
    "('element": {
      "type": "str",
      "description": "Element symbol (e.g., O, U, He)'),"
    }
  }
}
```
