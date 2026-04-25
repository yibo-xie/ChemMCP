# GetBondEnergy

> **Version:** 0.1.0 | **Function:** `get_bond_energy`

## Description
Query bond dissociation energies (BDE) in kJ/mol for common chemical bonds.

## Implementation Details
Uses a database of average bond dissociation energies at 298 K from standard references (CRC Handbook, NIST). Returns energy required to homolytically cleave a bond in the gas phase.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('bond` | str | N/A | Bond specification (e.g., C-C, C=C, C≡C, C-H, O-H, N≡N)'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `bond` | str | 'Bond specification' |
| `energy_kj_mol` | float | 'Bond dissociation energy in kJ/mol' |
| `energy_kcal_mol` | float | 'Bond dissociation energy in kcal/mol' |
| `note` | str | 'Interpretation note' |

## Examples
```python
examples = [
{'code_input': {'bond': 'C-C'}, 'text_input': {'bond': 'C-C'}, 'output': {'energy_kj_mol': 347, 'energy_kcal_mol': 83, 'note': 'Average C-C single bond', 'bond': 'C-C'}},
        {'code_input': {'bond': 'C=O'}, 'text_input': {'bond': 'C=O'}, 'output': {'energy_kj_mol': 799, 'energy_kcal_mol': 191, 'note': 'Carbonyl (formaldehyde/ketone)', 'bond': 'C=O'}},
        {'code_input': {'bond': 'N≡N'}, 'text_input': {'bond': 'N≡N'}, 'output': {'energy_kj_mol': 945, 'energy_kcal_mol': 226, 'note': 'N≡N triple bond (very strong)', 'bond': 'N≡N'}},
]
```

## Cherry Studio JSON Config
```json
{
  "name": "GetBondEnergy",
  "description": "Query bond dissociation energies (BDE) in kJ/mol for common chemical bonds.",
  "function": "get_bond_energy",
  "parameters": {
    "('bond": {
      "type": "str",
      "description": "Bond specification (e.g., C-C, C=C, C≡C, C-H, O-H, N≡N)'),"
    }
  }
}
```
