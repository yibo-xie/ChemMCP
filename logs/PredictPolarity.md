# PredictPolarity

> **Version:** 0.1.0 | **Function:** `predict_polarity`

## Description
Predict molecular polarity (polar/nonpolar) from SMILES string, including dipole moment direction and explanation.

## Implementation Details
Uses RDKit to analyze molecular structure and electronegativity differences between bonded atoms. Determines if bond dipoles cancel out due to symmetry. Returns polarity prediction, dipole moment estimate, and detailed reasoning.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('smiles` | str | N/A | SMILES string of the molecule'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `smiles` | str | 'Input SMILES' |
| `is_polar` | bool | 'Whether the molecule is polar' |
| `polarity` | str | '"polar" or "nonpolar"' |
| `dipole_analysis` | dict | 'Bond-by-bond dipole analysis' |
| `symmetry_analysis` | str | 'Symmetry and cancellation analysis' |
| `explanation` | str | 'Detailed explanation of polarity prediction' |

## Examples
```python
examples = [
{'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'smiles': 'CCO', 'is_polar': True, 'polarity': 'polar', 'dipole_analysis': [...
]
```

## Cherry Studio JSON Config
```json
{
  "name": "PredictPolarity",
  "description": "Predict molecular polarity (polar/nonpolar) from SMILES string, including dipole moment direction and explanation.",
  "function": "predict_polarity",
  "parameters": {
    "('smiles": {
      "type": "str",
      "description": "SMILES string of the molecule'),"
    }
  }
}
```
