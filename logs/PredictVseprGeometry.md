# PredictVseprGeometry

> **Version:** 0.1.0 | **Function:** `predict_vsepr_geometry`

## Description
Predict molecular geometry using VSEPR theory based on bonding pairs and lone pairs around a central atom.

## Implementation Details
Uses VSEPR (Valence Shell Electron Pair Repulsion) theory to predict electron geometry, molecular geometry, bond angles, and hybridization from steric number (bonding pairs + lone pairs).

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('molecule` | str | N/A | Molecular formula or SMILES string (e.g., H2O, NH3, SF6, C=O)'), |
| `('bonding_pairs` | int | N/A | Optional: number of bonding pairs around central atom'), |
| `('lone_pairs` | int | N/A | Optional: number of lone pairs on central atom'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `molecule` | str | 'Molecule identifier' |
| `steric_number` | int | 'Total number of electron domains (BP + LP |
| `electron_geometry` | str | 'Electron pair geometry' |
| `molecular_geometry` | str | 'Molecular geometry (atom positions only |
| `bond_angles` | str | 'Approximate bond angles' |
| `hybridization` | str | 'Predicted hybridization of central atom' |
| `description` | str | 'Detailed geometry description' |

## Examples
```python
examples = [
{'code_input': {'molecule': 'H2O', 'bonding_pairs': 2, 'lone_pairs': 2}, 'text_input': {'query': 'H2O'}, 'output': {'molecule': 'H2O', 'steric_number': 4, 'electron_geometry': 'tetrahedral', 'molecular_geometry': 'bent/angular', 'bond_angles': '104.5', 'hybridization': 'sp3', 'description': '...'}},
        {'code_input': {'molecule': 'SF6', 'bonding_pairs': 6, 'lone_pairs': 0}, 'text_input': {'query': 'SF6'}, 'output': {'molecule': 'SF6', 'steric_number': 6, 'electron_geometry': 'octahedral', 'molecular_geometry': 'octahedral', 'bond_angles': '90', 'hybridization': 'sp3d2', 'description': '...'}},
]
```

## Cherry Studio JSON Config
```json
{
  "name": "PredictVseprGeometry",
  "description": "Predict molecular geometry using VSEPR theory based on bonding pairs and lone pairs around a central atom.",
  "function": "predict_vsepr_geometry",
  "parameters": {
    "('molecule": {
      "type": "str",
      "description": "Molecular formula or SMILES string (e.g., H2O, NH3, SF6, C=O)'),"
    },
    "('bonding_pairs": {
      "type": "int",
      "description": "Optional: number of bonding pairs around central atom'),"
    },
    "('lone_pairs": {
      "type": "int",
      "description": "Optional: number of lone pairs on central atom'),"
    }
  }
}
```
