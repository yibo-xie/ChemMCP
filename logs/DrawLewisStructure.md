# DrawLewisStructure

> **Version:** 0.1.0 | **Function:** `draw_lewis_structure`

## Description
Draw the Lewis structure of a molecule showing valence electrons, bonding pairs, and lone pairs.

## Implementation Details
Uses a rule-based algorithm with predefined Lewis structures for common molecules. Shows electron dot notation, bond orders, lone pairs, octet status, and molecular geometry description.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('formula` | str | N/A | Molecular formula (e.g., H2O, CO2, NH3, CH4, SO2, PCl5, SF6)'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `formula` | str | 'Molecular formula' |
| `lewis_diagram` | str | 'ASCII/text Lewis structure diagram' |
| `total_valence_electrons` | int | 'Total count of valence electrons' |
| `bond_info` | list | 'List of bonds between atoms' |
| `lone_pair_info` | dict | 'Lone pair distribution per atom' |
| `geometry_description` | str | 'Molecular geometry description' |
| `octet_status` | str | 'Whether octet rule is satisfied' |

## Examples
```python
examples = [
{'code_input': {'formula': 'H2O'}, 'text_input': {'formula': 'H2O'}, 'output': {'formula': 'H2O', 'lewis_diagram': '...', 'total_valence_electrons': 8, 'bond_info': '...', 'lone_pair_info': '...', 'geometry_description': 'bent (104.5)', 'octet_status': '...'}},
        {'code_input': {'formula': 'CO2'}, 'text_input': {'formula': 'CO2'}, 'output': {'formula': 'CO2', 'lewis_diagram': '...', 'total_valence_electrons': 16, 'bond_info': '...', 'lone_pair_info': '...', 'geometry_description': 'linear (180)', 'octet_status': '...'}},
        {'code_input': {'formula': 'SF6'}, 'text_input': {'formula': 'SF6'}, 'output': {'formula': 'SF6', 'lewis_diagram': '...', 'total_valence_electrons': 48, 'bond_info': '...', 'lone_pair_info': '...', 'geometry_description': 'octahedral (90)', 'octet_status': '...'}},
]
```

## Cherry Studio JSON Config
```json
{
  "name": "DrawLewisStructure",
  "description": "Draw the Lewis structure of a molecule showing valence electrons, bonding pairs, and lone pairs.",
  "function": "draw_lewis_structure",
  "parameters": {
    "('formula": {
      "type": "str",
      "description": "Molecular formula (e.g., H2O, CO2, NH3, CH4, SO2, PCl5, SF6)'),"
    }
  }
}
```
