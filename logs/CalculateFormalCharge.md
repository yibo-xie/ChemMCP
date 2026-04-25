# CalculateFormalCharge

> **Version:** 0.1.0 | **Function:** `calculate_formal_charge`

## Description
Calculate formal charge for each atom in a molecule from its SMILES string.

## Implementation Details
Uses RDKit to parse the molecule and calculates formal charge as FC = V - N - B/2, where V is valence electrons of neutral atom, N is number of non-bonding electrons, and B is number of bonding electrons. Returns per-atom charges and total molecular charge.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('smiles` | str | N/A | SMILES string of the molecule'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `atom_charges` | list | 'Formal charge for each atom with atom index, symbol, and charge' |
| `total_charge` | int | 'Total molecular charge (sum of all atom charges |
| `smiles` | str | 'Input SMILES' |

## Examples
```python
examples = [
{'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'atom_charges': [...
]
```

## Cherry Studio JSON Config
```json
{
  "name": "CalculateFormalCharge",
  "description": "Calculate formal charge for each atom in a molecule from its SMILES string.",
  "function": "calculate_formal_charge",
  "parameters": {
    "('smiles": {
      "type": "str",
      "description": "SMILES string of the molecule'),"
    }
  }
}
```
