# AnalyzeMolecularOrbital

> **Version:** 0.1.0 | **Function:** `analyze_molecular_orbital`

## Description
Analyze molecular orbital diagram for simple diatomic molecules: bond order, magnetic property, MO configuration.

## Implementation Details
Implements molecular orbital theory for period 2 homonuclear diatomic molecules (H2 through Ne2) and common heteronuclear diatomics (CO, NO, HF, CN). Returns electron configuration, bond order, magnetic behavior, and stability analysis.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('molecule` | str | N/A | Diatomic molecule formula (e.g., N2, O2, CO, NO, HF)'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `molecule` | str | 'Molecule identifier' |
| `electron_configuration` | str | 'MO electron configuration' |
| `bond_order` | float | 'Bond order (number of chemical bonds |
| `magnetic_property` | str | 'Diamagnetic or paramagnetic' |
| `stability` | str | 'Stability assessment' |
| `description` | str | 'Detailed explanation of MO analysis' |

## Examples
```python
examples = [
{'code_input': {'molecule': 'N2'}, 'text_input': {'molecule': 'N2'}, 'output': {'molecule': 'N2', 'electron_configuration': '...', 'bond_order': 3, 'magnetic_property': 'diamagnetic', 'stability': 'very stable', 'description': '...'}},
        {'code_input': {'molecule': 'O2'}, 'text_input': {'molecule': 'O2'}, 'output': {'molecule': 'O2', 'electron_configuration': '...', 'bond_order': 2, 'magnetic_property': 'paramagnetic', 'stability': 'stable', 'description': '...'}},
        {'code_input': {'molecule': 'CO'}, 'text_input': {'molecule': 'CO'}, 'output': {'molecule': 'CO', 'electron_configuration': '...', 'bond_order': 3, 'magnetic_property': 'diamagnetic', 'stability': 'stable', 'description': '...'}},
]
```

## Cherry Studio JSON Config
```json
{
  "name": "AnalyzeMolecularOrbital",
  "description": "Analyze molecular orbital diagram for simple diatomic molecules: bond order, magnetic property, MO configuration.",
  "function": "analyze_molecular_orbital",
  "parameters": {
    "('molecule": {
      "type": "str",
      "description": "Diatomic molecule formula (e.g., N2, O2, CO, NO, HF)'),"
    }
  }
}
```
