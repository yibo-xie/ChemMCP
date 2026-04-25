# ElementAbundance

> **Version:** 0.1.0 | **Function:** `element_abundance`

## Description
Query element abundance in Earth's crust, oceans, and universe.

## Implementation Details
Returns abundance data from geochemical surveys and astrophysical measurements. Crust abundance in ppm (parts per million by mass), ocean concentration in mg/L, cosmic abundance relative to silicon = 10^6.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('element` | str | N/A | Element symbol (e.g., O, Fe, Au)'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `element` | str | 'Element symbol' |
| `crust_abundance_ppm` | float | 'Abundance in Earth\'s crust (ppm |
| `ocean_concentration_mg_L` | float | 'Concentration in ocean water (mg/L |
| `cosmic_abundance` | float | 'Cosmic abundance (Si=1e6 |
| `rank_note` | str | 'Ranking context' |

## Examples
```python
examples = [
{'code_input': {'element': 'Fe'}, 'text_input': {'element': 'Fe'}, 'output': {'element': 'Fe', 'crust_abundance_ppm': 56300, 'ocean_concentration_mg_L': None, 'cosmic_abundance': None, 'rank_note': '...'}},
        {'code_input': {'element': 'Au'}, 'text_input': {'element': 'Au'}, 'output': {'element': 'Au', 'crust_abundance_ppm': 0.004, 'ocean_concentration_mg_L': None, 'cosmic_abundance': None, 'rank_note': '...'}},
]
```

## Cherry Studio JSON Config
```json
{
  "name": "ElementAbundance",
  "description": "Query element abundance in Earth's crust, oceans, and universe.",
  "function": "element_abundance",
  "parameters": {
    "('element": {
      "type": "str",
      "description": "Element symbol (e.g., O, Fe, Au)'),"
    }
  }
}
```
